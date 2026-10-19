#!/usr/bin/perl
use strict;
use warnings;
use File::Basename;

# MATLAB to C++/Doxygen Converter
# Converts MATLAB files to C++ with proper Doxygen documentation
# Following Richard K. Johnson's MATLAB Style Guidelines

if (@ARGV != 1) {
    die "Usage: $0 <matlab_file.m>\n";
}

my $filename = $ARGV[0];
my $output = "";

# Extract package information from file path
my ($base_name, $dir_path, $suffix) = fileparse($filename, qr/\.[^.]*/);
my $package_info = extract_package_info($dir_path);

# Parse the MATLAB file
my $file_content = parse_matlab_file($filename);

# Generate C++ output with Doxygen documentation
$output = generate_cpp_output($file_content, $base_name, $package_info);

print $output;

#==============================================================================
# SUBROUTINES
#==============================================================================

sub parse_matlab_file {
    my ($filename) = @_;
    
    open(my $fh, '<', $filename) or die "Cannot open file $filename: $!";
    my @lines = <$fh>;
    close($fh);
    
    my %file_data = (
        'docstring' => [],
        'function_signature' => '',
        'class_signature' => '',
        'is_class' => 0,
        'raw_lines' => \@lines,
        'processed_content' => ''
    );
    
    my $current_section = 'none';
    my $docstring_started = 0;
    my $in_docstring = 0;
    
    # First pass: extract docstring and determine type
    for my $i (0..$#lines) {
        my $line = $lines[$i];
        chomp $line;
        
        # Skip empty lines at the beginning
        next if $line =~ /^\s*$/ && !$docstring_started;
        
        # Detect function declaration
        if ($line =~ /^\s*function\s+(.*?)$/) {
            $file_data{'function_signature'} = $1;
            $file_data{'is_class'} = 0;
            $current_section = 'function';
            
            # Look for docstring starting from next line
            for my $j ($i+1..$#lines) {
                my $next_line = $lines[$j];
                chomp $next_line;
                
                if ($next_line =~ /^\s*%(?!<)(.*)$/) {
                    push @{$file_data{'docstring'}}, $1;
                    $in_docstring = 1;
                    $docstring_started = 1;
                } elsif ($in_docstring && $next_line =~ /^\s*$/) {
                    # Allow empty lines in docstring
                    push @{$file_data{'docstring'}}, '';
                } elsif ($in_docstring) {
                    # End of docstring
                    last;
                } elsif ($next_line =~ /^\s*$/) {
                    # Skip empty lines before docstring
                    next;
                } else {
                    # No docstring found
                    last;
                }
            }
            last; # We found the function, stop parsing
        }
        
        # Detect class declaration
        elsif ($line =~ /^\s*classdef\s*(.*?)$/) {
            $file_data{'class_signature'} = $1;
            $file_data{'is_class'} = 1;
            $current_section = 'class';
            
            # Look for docstring starting from next line
            for my $j ($i+1..$#lines) {
                my $next_line = $lines[$j];
                chomp $next_line;
                
                if ($next_line =~ /^\s*%(.*)$/) {
                    push @{$file_data{'docstring'}}, $1;
                    $in_docstring = 1;
                    $docstring_started = 1;
                } elsif ($in_docstring && $next_line =~ /^\s*$/) {
                    # Allow empty lines in docstring
                    push @{$file_data{'docstring'}}, '';
                } elsif ($in_docstring) {
                    # End of docstring
                    last;
                } elsif ($next_line =~ /^\s*$/) {
                    # Skip empty lines before docstring
                    next;
                } else {
                    # No docstring found
                    last;
                }
            }
            last; # We found the class, stop parsing for docstring
        }
    }
    
    return \%file_data;
}

sub prepare_content_without_docstring {
    my ($file_data) = @_;
    
    my @lines = @{$file_data->{'raw_lines'}};
    my @content_lines = ();
    my $found_classdef = 0;
    my $skip_docstring = 0;
    my $docstring_ended = 0;
    
    for my $i (0..$#lines) {
        my $line = $lines[$i];
        chomp $line;
        
        # Find classdef line
        if ($line =~ /^\s*classdef\s/) {
            $found_classdef = 1;
            push @content_lines, $line;
            $skip_docstring = 1; # Start skipping docstring after classdef
            next;
        }
        
        # If we found classdef and are in docstring skipping mode
        if ($found_classdef && $skip_docstring && !$docstring_ended) {
            # Skip empty lines immediately after classdef
            if ($line =~ /^\s*$/) {
                next;
            }
            # Skip comment lines (docstring)
            elsif ($line =~ /^\s*%/) {
                next;
            }
            # First non-comment, non-empty line ends docstring
            else {
                $docstring_ended = 1;
                $skip_docstring = 0;
                push @content_lines, $line;
            }
        }
        else {
            # Include all other lines
            push @content_lines, $line;
        }
    }
    
    return join("\n", @content_lines);
}

# Helper function to check if a method is a constructor
sub is_constructor {
    my ($method_name, $class_name) = @_;
    return defined($class_name) && $method_name eq $class_name;
}

# Helper function to check if a method is a built-in get/set method
sub is_builtin_getter_setter {
    my ($method_name) = @_;
    # Match patterns like get.propertyName or set.propertyName
    return $method_name =~ /^(get|set)\./;
}

sub convert_class_content {
    my ($content_text, $class_signature) = @_;
    
    my @lines = split /\n/, $content_text;
    
    my $output = "";
    my $declTypeDef = "";
    my $inClass = 0;
    my $inAbstractMethodBlock = 0;
    my $listeProperties = 0;
    my $listeEnumeration = 0;
    my $listeEvents = 0;
    my $methodAttribute = "";
    my $typeProperties = "";
    my $className = "";
    my $inMethodDocstring = 0;
    my @method_docstring = ();
    my $current_function_name = "";
    my $classdef_processed = 0;
    my $pending_method = "";
    
    # Extract class name from signature for constructor detection
    if ($class_signature =~ /^\s*(\w+)/) {
        $className = $1;
    }
    
    for my $line (@lines) {
        chomp $line;
        
        # Handle comment conversion (including method docstrings)
        if ($line =~ /(^\s*)(%>)(.*)/) {
            $output .= "$1///$3\n";
            next;
        }
        
        # Handle regular comments in methods (potential docstrings)
        if ($line =~ /^\s*%(?!<)(.*)$/ && $inClass) {
            my $comment_content = $1;
            if ($inMethodDocstring) {
                push @method_docstring, $comment_content;
            } else {
                # Start of a method docstring
                $inMethodDocstring = 1;
                @method_docstring = ($comment_content);
            }
            next;
        } elsif ($inMethodDocstring && $line !~ /^\s*%/ && $line !~ /^\s*$/) {
            # End of method docstring, process it
            my $method_doc = generate_method_docstring(\@method_docstring, $current_function_name, $className);

            $output .= $method_doc;
            if ($pending_method){
                $output .= $pending_method;
                $pending_method = "";
            }
            $inMethodDocstring = 0;
            @method_docstring = ();
        }
        
        # Handle end of properties block
        if ($listeProperties == 1 && $line =~ /^\s*\bend\b\s*/) {
            $listeProperties = 0;
        }
        
        # Handle end of abstract method block
        if ($inAbstractMethodBlock == 1 && $line =~ /^\s*\bend\b\s*/) {
            $inAbstractMethodBlock = 0;
        }
        
        # Handle properties within properties block (preserve comments)
        if ($listeProperties == 1 && $line =~ /^\s*([\w\d]*)\s*([^%]*)?(.*)/) {
            my $propertyName = $1;
            my $propertyDeclaration = $2 || '';
            my $propertyComment = $3 || '';
            
            # Clean up the property declaration (remove extra whitespace)
            $propertyDeclaration =~ s/^\s+|\s+$//g;
            
            if ($propertyName !~ /^$/) {
                my $properties;
                if ($typeProperties =~ /Constant/) {
                    $properties = $propertyName . " " . $propertyDeclaration . ";";
                } else {
                    $properties = $propertyName . " " . $propertyDeclaration . ";";
                }
                
                # Clean up extra spaces in the property declaration
                $properties =~ s/\s+/ /g;
                $properties =~ s/\s+;/;/g;
                
                # Handle property comments - only include if NOT following "PROPERTYNAME - DESCRIPTION" pattern
                if ($propertyComment =~ /^\s*%\s*(.*)$/) {
                    my $comment = $1;
                    $comment =~ s/^\s+|\s+$//g;
                    
                    # Only include comment if it does NOT start with property name followed by dash
                    if (!($comment =~ /^$propertyName\s*-\s*(.+)$/)) {
                        # Use ///< for same-line property comments
                        $properties .= " ///< $comment";
                    }
                    # If it follows PROPERTYNAME - DESCRIPTION pattern, skip the comment entirely
                }
                
                $output .= $typeProperties . "Property " . $properties . "\n";
            }
        }
        
        # Handle end of enumeration block
        if ($listeEnumeration == 1 && $line =~ /^\s*\bend\b\s*/) {
            $listeEnumeration = 0;
            $output .= "};\n";
        }
        
        # Handle end of events block
        if ($listeEvents == 1 && $line =~ /^\s*\bend\b\s*/) {
            $listeEvents = 0;
            $output .= "};\n";
        }
        
        # Handle events within events block
        if ($listeEvents == 1 && $line =~ /^\s*([\w\d]*)\s*(.*)/) {
            my $name_event = $1;
            my $comment = $2 || '';
            if ($name_event !~ /^$/) {
                my $event = $name_event . ",";
                if ($comment =~ /^\s*%(.*)$/) {
                    my $event_comment = $1;
                    $event_comment =~ s/^\s+|\s+$//g;
                    $event .= " /// $event_comment";
                }
                $output .= $event . "\n";
            }
        }
        
        # Handle enumeration values within enumeration block
        if ($listeEnumeration == 1 && $line =~ /^\s*([\w\d]*)\s*(\(.*\))?\s*(.*)/) {
            my $name_enum = $1;
            my $val_enum = $2 || '';
            my $comment = $3 || '';
            
            if ($name_enum !~ /^$/) {
                my $enum;
                if ($val_enum !~ /^$/) {
                    $enum = "$name_enum=$val_enum,";
                } else {
                    $enum = "$name_enum,";
                }
                
                if ($comment =~ /^\s*%(.*)$/) {
                    my $enum_comment = $1;
                    $enum_comment =~ s/^\s+|\s+$//g;
                    $enum .= " /// $enum_comment";
                }
                
                $output .= $enum . "\n";
            }
        }
        
        # Handle function declarations
        if ($line =~ /(^\s*function)\s*([\] \w\d,_\[]+=)?\s*([.\w\d_-]*)\s*\(?([\w\d\s,~]*)\)?(.*)/) {
            my $functionKeyWord = $1;
            my $returns = $2 || '';
            my $functionName = $3;
            my $arguments = $4 || '';
            
            $current_function_name = $functionName;
            
            # Skip built-in get/set methods entirely
            if (is_builtin_getter_setter($functionName)) {
                next;
            }
            
            if ($inClass == 0) {
                $output = $declTypeDef . $output;
                $declTypeDef = "";
            }
            
            # Convert arguments with proper handling of varargin
            my @arg_parts = split /\s*,\s*/, $arguments;
            my @converted_args = ();
            
            for my $arg (@arg_parts) {
                $arg =~ s/^\s+|\s+$//g; # trim
                if ($arg eq '') {
                    next;
                } elsif ($arg eq '~') {
                    push @converted_args, 'ignoredArg';
                } elsif ($arg eq 'varargin') {
                    push @converted_args, 'VariableInputType varargin';
                } else {
                    push @converted_args, "InputType $arg";
                }
            }
            
            $arguments = join(', ', @converted_args);
            
            my $ligne = "$methodAttribute $functionKeyWord $functionName($arguments);\n";
            $pending_method .= $ligne;
        }
        # Handle abstract method signatures
        elsif ($line =~ /^\s*([\] \w\d,_\[]+=)?\s*([.\w\d_-]+)\s*\(?([\w\d\s,~]*)\)?(.*)/ && $inAbstractMethodBlock == 1) {
            my $functionName = $2;
            my $arguments = $3 || '';
            
            $current_function_name = $functionName;
            
            # Skip built-in get/set methods entirely
            if (is_builtin_getter_setter($functionName)) {
                next;
            }
            
            # Convert arguments similar to above
            my @arg_parts = split /\s*,\s*/, $arguments;
            my @converted_args = ();
            
            for my $arg (@arg_parts) {
                $arg =~ s/^\s+|\s+$//g;
                if ($arg eq '') {
                    next;
                } elsif ($arg eq '~') {
                    push @converted_args, 'ignoredArg';
                } elsif ($arg eq 'varargin') {
                    push @converted_args, 'VariableInputType varargin';
                } else {
                    push @converted_args, "InputType $arg";
                }
            }
            
            $arguments = join(', ', @converted_args);
            
            my $ligne = "$methodAttribute function $functionName($arguments);\n";
            $output .= $ligne;
        }
        
        # Handle class definition - but only process variables, don't output the definition again
        if ($line =~ /(^\s*classdef)\s*(\s*\([\{\}\?\w,=\s]+\s*\))?\s*([\w\d_]+)\s*<?\s*([\s\w\d._&]+)?(.*)/) {
            if (!$classdef_processed) {
                $className = $3;
                my $classInheritance = $4 || '';
                my $classAttributes = $2 || '';
                
                # Generate class definition from passed signature instead of parsing again
                my $classDef = generate_class_definition_from_signature($class_signature);
                $output .= $classDef;
                $inClass = 1;
                $classdef_processed = 1;
            }
            next; # Skip this line since we've already processed the class definition
        }
        
        # Handle properties block
        if ($line =~ /(^\s*properties)\s*(\s*\([\w,=\s]+\s*\))?(.*)/) {
            $listeProperties = 1;
            my $propertiesAttributes = $2 || '';
            $typeProperties = "public:\n";
            
            if (lc($propertiesAttributes) =~ /(access\s*=\s*private)/) {
                $typeProperties = "private:\n";
            } elsif (lc($propertiesAttributes) =~ /(access\s*=\s*public)/) {
                $typeProperties = "public:\n";
            } elsif (lc($propertiesAttributes) =~ /(access\s*=\s*protected)/) {
                $typeProperties = "protected:\n";
            }
            
            if (lc($propertiesAttributes) =~ /(constant(\s*=\s*true\s*)?)/ && 
                !(lc($propertiesAttributes) =~ /(constant\s*=\s*false)/ || lc($propertiesAttributes) =~ /(~constant)/)) {
                $typeProperties = $typeProperties . " Constant ";
            }
        }
        
        # Handle enumeration block
        if ($line =~ /(^\s*enumeration)\s*(.*)/) {
            $listeEnumeration = 1;
            $output .= "public:\nenum $className {\n";
        }
        
        # Handle events block
        if ($line =~ /(^\s*events)\s*(.*)/) {
            $listeEvents = 1;
            $output .= "public:\nenum Events {\n";
        }
        
        # Handle methods block
        if ($line =~ /(^\s*methods)\s*(\s*\([\w,=\s]+\s*\))?(.*)/) {
            $methodAttribute = "public:\n";
            my $methodsAttributes = $2 || '';
            
            if (lc($methodsAttributes) =~ /(access\s*=\s*private)/) {
                $methodAttribute = "private:\n";
            } elsif (lc($methodsAttributes) =~ /(access\s*=\s*protected)/) {
                $methodAttribute = "protected:\n";
            } elsif (lc($methodsAttributes) =~ /(access\s*=\s*public)/) {
                $methodAttribute = "public:\n";
            }
            
            if (lc($methodsAttributes) =~ /(abstract(\s*=\s*true\s*)?)/) {
                $inAbstractMethodBlock = 1;
                $methodAttribute = $methodAttribute . " virtual ";
            }
            
            if (lc($methodsAttributes) =~ /(static(\s*=\s*true\s*)?)/ && 
                !(lc($methodsAttributes) =~ /(static\s*=\s*false)/ || lc($methodsAttributes) =~ /(~static)/)) {
                $methodAttribute = $methodAttribute . " static";
            }
        }
        
        # Only add newline if we haven't skipped the classdef line
        if (!($line =~ /^\s*classdef/ && $classdef_processed)) {
            $output .= "\n";
        }
    }
    
    # Close class if it was opened
    if ($inClass) {
        $output .= "};\n";
    }
    
    return $output;
}

sub generate_class_definition_from_signature {
    my ($class_signature) = @_;
    
    # Parse the class signature to generate proper C++ class definition
    # Remove abstract attribute from class definition
    if ($class_signature =~ /^\s*(\([^)]*\))?\s*(\w+)(?:\s*<\s*(.+))?\s*$/) {
        my $attributes = $1 || '';
        my $class_name = $2;
        my $inheritance = $3;
        
        my $class_def = "class $class_name";
        if ($inheritance) {
            $inheritance =~ s/&/, public /g;
            $class_def .= " : public $inheritance";
        }
        $class_def .= " {\npublic:\n";
        return $class_def;
    }
    
    # Fallback for simple class name without attributes
    if ($class_signature =~ /^\s*(\w+)(?:\s*<\s*(.+))?\s*$/) {
        my $class_name = $1;
        my $inheritance = $2;
        
        my $class_def = "class $class_name";
        if ($inheritance) {
            $inheritance =~ s/&/, public /g;
            $class_def .= " : public $inheritance";
        }
        $class_def .= " {\npublic:\n";
        return $class_def;
    }
    
    # Fallback
    return "class UnknownClass {\npublic:\n";
}

sub generate_method_docstring {
    my ($docstring_lines, $method_name, $class_name) = @_;
    
    return "" if !@$docstring_lines;
    
    # Check if this is a built-in get/set method - if so, skip entirely
    if (is_builtin_getter_setter($method_name)) {
        return "";
    }
    
    my $output = "";
    
    # Parse method docstring using the updated parse_docstring function
    my $parsed_doc = parse_docstring($docstring_lines);
    
    # Add brief description
    if ($parsed_doc->{'brief'}) {
        $output .= "/// \@brief " . $parsed_doc->{'brief'} . "\n";
    }
    
    # Add detailed description
    if ($parsed_doc->{'detailed'}) {
        $output .= "///\n";
        my @detail_lines = split /\n/, $parsed_doc->{'detailed'};
        
        for my $line (@detail_lines) {
            $line =~ s/^\s+|\s+$//g; # trim
            if ($line ne '') {
                $output .= "/// $line\n";
            } else {
                $output .= "///\n";
            }
        }
    }
    
    # Parameters (inputs)
    if (@{$parsed_doc->{'inputs'}}) {
        $output .= "///\n";
        for my $input (@{$parsed_doc->{'inputs'}}) {
            if ($input =~ /^(\w+)\s*-\s*(.+)$/) {
                $output .= "/// \@param $1 - $2\n";
            } else {
                $output .= "/// \@param $input\n";
            }
        }
    }
    
    # Return values (outputs) - skip for constructors
    my $is_constructor_method = is_constructor($method_name, $class_name);
    if (@{$parsed_doc->{'outputs'}} && !$is_constructor_method) {
        $output .= "///\n";
        for my $output_item (@{$parsed_doc->{'outputs'}}) {
            if ($output_item =~ /^(\w+)\s*-\s*(.+)$/) {
                $output .= "/// \@retval $1 - $2\n";
            } else {
                $output .= "/// \@retval $output_item\n";
            }
        }
    }
    
    # Examples
    if (@{$parsed_doc->{'examples'}}) {
        $output .= "///\n";
        $output .= "/// \@code\n";
        for my $example (@{$parsed_doc->{'examples'}}) {
            $output .= "/// $example\n";
        }
        $output .= "/// \@endcode\n";
    }
    
    # Notes
    if ($parsed_doc->{'notes'}) {
        $output .= "///\n";
        for my $line (split /\n/, $parsed_doc->{'notes'}) {
            $output .= "/// \@note $line\n";
        }
    }
    
    $output .= "///\n";
    return $output;
}

sub extract_package_info {
    my ($dir_path) = @_;
    
    # Extract package hierarchy from directory path
    # Example: /path/to/+core/+linalg/ becomes core::linalg
    my @package_parts = ();
    
    while ($dir_path =~ /\+(\w+)/g) {
        push @package_parts, $1;
    }
    
    return {
        'namespace' => join('::', @package_parts),
        'parts' => \@package_parts
    };
}

sub parse_docstring {
    my ($docstring_lines) = @_;
    
    my %parsed = (
        'brief' => '',
        'detailed' => '',
        'inputs' => [],
        'outputs' => [],
        'examples' => [],
        'notes' => '',
        'see_also' => []
    );
    
    return \%parsed if !@$docstring_lines;
    
    # Parse first line for FUNCTIONNAME SUMMARY pattern and handle multi-line summary
    my $first_line = $docstring_lines->[0] || '';
    $first_line =~ s/^\s+|\s+$//g; # trim whitespace
    
    my $function_name = '';
    my $brief_content = '';
    my $brief_end_index = 0;
    
    if ($first_line =~ /^(\w+)\s+(.+)$/) {
        $function_name = $1;
        $brief_content = $2;
        $parsed{'function_name'} = $function_name;
        
        # Look for continuation of brief summary in subsequent lines
        for my $i (1..$#{$docstring_lines}) {
            my $line = $docstring_lines->[$i];
            $line =~ s/^\s+|\s+$//g; # trim whitespace
            
            # Stop if we hit an empty line (end of brief)
            if ($line eq '') {
                $brief_end_index = $i;
                last;
            }
            
            # Stop if we hit a section header
            if ($line =~ /^(Inputs?|Outputs?|Examples?|Notes?|See\s+Also):\s*$/i) {
                $brief_end_index = $i - 1;
                last;
            }
            
            # Otherwise, this line continues the brief summary
            $brief_content .= " " . $line;
            $brief_end_index = $i;
        }
        
        $parsed{'brief'} = $brief_content;
    } else {
        $parsed{'brief'} = $first_line;
        $brief_end_index = 0;
    }
    
    my $current_section = 'detailed';
    my $detailed_started = 0;
    
    # Start parsing from after the brief section
    for my $i (($brief_end_index + 1)..$#{$docstring_lines}) {
        my $line = $docstring_lines->[$i];
        $line =~ s/^\s+|\s+$//g; # trim whitespace
        
        # Skip empty lines
        if ($line eq '') {
            # If we're in detailed section and haven't started yet, this might be the separator
            if ($current_section eq 'detailed' && !$detailed_started) {
                next;
            }
            # If we're already in detailed section, preserve empty lines
            elsif ($current_section eq 'detailed') {
                $parsed{'detailed'} .= "\n";
            }
            next;
        }
        
        # Detect section headers (syntax section removed)
        if ($line =~ /^Inputs?:\s*$/i) {
            $current_section = 'inputs';
            next;
        } elsif ($line =~ /^Outputs?:\s*$/i) {
            $current_section = 'outputs';
            next;
        } elsif ($line =~ /^Examples?:\s*$/i) {
            $current_section = 'examples';
            next;
        } elsif ($line =~ /^Notes?:\s*$/i) {
            $current_section = 'notes';
            next;
        } elsif ($line =~ /^See\s+Also:\s*$/i) {
            $current_section = 'see_also';
            next;
        }
        
        # Process content based on current section
        if ($current_section eq 'detailed') {
            if (!$detailed_started) {
                $parsed{'detailed'} = $line;
                $detailed_started = 1;
            } else {
                $parsed{'detailed'} .= "\n" . $line;
            }
        } elsif ($current_section eq 'inputs') {
            push @{$parsed{'inputs'}}, $line;
        } elsif ($current_section eq 'outputs') {
            push @{$parsed{'outputs'}}, $line;
        } elsif ($current_section eq 'examples') {
            push @{$parsed{'examples'}}, $line;
        } elsif ($current_section eq 'notes') {
            if ($parsed{'notes'} eq '') {
                $parsed{'notes'} = $line;
            } else {
                $parsed{'notes'} .= "\n" . $line;
            }
        } elsif ($current_section eq 'see_also') {
            # Split comma-separated references
            my @refs = split /,\s*/, $line;
            push @{$parsed{'see_also'}}, @refs;
        }
    }
    
    # Clean up trailing newlines
    $parsed{'detailed'} =~ s/\n+$//;
    $parsed{'notes'} =~ s/\n+$//;
    
    return \%parsed;
}

sub generate_file_header_comment {
    my ($parsed_doc, $function_name, $package_info) = @_;
    
    my $header = "";
    
    # Add file documentation at the top - only @file and @namespace
    if ($function_name) {
        $header .= "/// \@file $function_name.m\n";
        if ($package_info->{'namespace'}) {
            $header .= "/// \@namespace $package_info->{'namespace'}\n";
        }
        $header .= "///\n";
    }
    
    return $header;
}

sub generate_function_comment {
    my ($parsed_doc, $package_info) = @_;
    
    my $doxygen = "";
    
    # Add brief from first line summary
    if ($parsed_doc->{'brief'}) {
        $doxygen .= "/// \@brief " . $parsed_doc->{'brief'} . "\n";
    }
    
    # Add detailed description from second paragraph
    if ($parsed_doc->{'detailed'}) {
        $doxygen .= "///\n";
        my @detail_lines = split /\n/, $parsed_doc->{'detailed'};
        
        for my $line (@detail_lines) {
            $line =~ s/^\s+|\s+$//g; # trim
            $doxygen .= "/// $line\n";
        }
    }
    
    # Parameters (inputs)
    if (@{$parsed_doc->{'inputs'}}) {
        $doxygen .= "///\n";
        for my $input (@{$parsed_doc->{'inputs'}}) {
            if ($input =~ /^(\w+)\s*-\s*(.+)$/) {
                $doxygen .= "/// \@param $1 - $2\n";
            } else {
                $doxygen .= "/// \@param $input\n";
            }
        }
    }
    
    # Return values (outputs) - using @retval instead of @return
    if (@{$parsed_doc->{'outputs'}}) {
        $doxygen .= "///\n";
        for my $output (@{$parsed_doc->{'outputs'}}) {
            if ($output =~ /^(\w+)\s*-\s*(.+)$/) {
                $doxygen .= "/// \@retval $1 - $2\n";
            } else {
                $doxygen .= "/// \@retval $output\n";
            }
        }
    }
    
    # Examples
    if (@{$parsed_doc->{'examples'}}) {
        $doxygen .= "///\n";
        $doxygen .= "/// \@code\n";
        for my $example (@{$parsed_doc->{'examples'}}) {
            $doxygen .= "/// $example\n";
        }
        $doxygen .= "/// \@endcode\n";
    }
    
    # Notes
    if ($parsed_doc->{'notes'}) {
        $doxygen .= "///\n";
        for my $line (split /\n/, $parsed_doc->{'notes'}) {
            $doxygen .= "/// \@note $line\n";
        }
    }
    
    # See also references with cross-reference conversion
    if (@{$parsed_doc->{'see_also'}}) {
        $doxygen .= "///\n";
        for my $ref (@{$parsed_doc->{'see_also'}}) {
            $ref =~ s/^\s+|\s+$//g; # trim
            my $converted_ref = convert_matlab_reference_to_cpp($ref, $package_info);
            $doxygen .= "/// \@see $converted_ref\n";
        }
    }
    
    return $doxygen;
}

sub convert_matlab_reference_to_cpp {
    my ($matlab_ref, $package_info) = @_;
    
    # Keep MATLAB dot notation in the output
    # Handle MATLAB package references like "core.linalg.solver" 
    # Keep as "core.linalg.solver" (no conversion to ::)
    if ($matlab_ref =~ /^([a-zA-Z_]\w*(?:\.[a-zA-Z_]\w*)*)\.([a-zA-Z_]\w*)$/) {
        my $package_path = $1;
        my $function_name = $2;
        
        return "$package_path.$function_name";
    }
    # Handle simple function references within same package
    elsif ($matlab_ref =~ /^[a-zA-Z_]\w*$/) {
        if ($package_info->{'namespace'}) {
            return "$package_info->{'namespace'}.$matlab_ref";
        } else {
            return $matlab_ref;
        }
    }
    # Handle class method references like "ClassName.methodName"
    elsif ($matlab_ref =~ /^([a-zA-Z_]\w*)\.([a-zA-Z_]\w*)$/) {
        my $class_name = $1;
        my $method_name = $2;
        
        if ($package_info->{'namespace'}) {
            return "$package_info->{'namespace'}.$class_name.$method_name";
        } else {
            return "$class_name.$method_name";
        }
    }
    
    # Fallback: return as-is if no pattern matches
    return $matlab_ref;
}

sub convert_matlab_signature_to_cpp {
    my ($signature, $is_class) = @_;
    
    if ($is_class) {
        # Handle class signature - following original m2cpp.pl pattern
        if ($signature =~ /^\s*(\w+)(?:\s*<\s*(.+))?\s*$/) {
            my $class_name = $1;
            my $inheritance = $2;
            
            my $cpp_sig = "class $class_name";
            if ($inheritance) {
                $inheritance =~ s/&/, public /g;
                $cpp_sig .= " : public $inheritance";
            }
            $cpp_sig .= " {\npublic:\n";
            return $cpp_sig;
        }
    } else {
        # Handle function signature - following original m2cpp.pl pattern
        if ($signature =~ /^(?:([\] \w\d,_\[]+=)\s*)?(\w+)\s*\(([\w\d\s,~]*)\)\s*$/) {
            my $returns = $1 || '';
            my $func_name = $2;
            my $arguments = $3 || '';
            
            # Handle special cases first: varargin and varargout
            $arguments =~ s/\bvarargin\b/VariableInputType varargin/g;
            
            # Convert other arguments with InputType, but skip already converted ones
            my @arg_parts = split /\s*,\s*/, $arguments;
            my @converted_args = ();
            
            for my $arg (@arg_parts) {
                $arg =~ s/^\s+|\s+$//g; # trim
                if ($arg eq '') {
                    # Skip empty arguments
                    next;
                } elsif ($arg eq '~') {
                    push @converted_args, 'ignoredArg';
                } elsif ($arg =~ /^VariableInputType/) {
                    # Already converted varargin
                    push @converted_args, $arg;
                } else {
                    # Regular argument
                    push @converted_args, "InputType $arg";
                }
            }
            
            $arguments = join(', ', @converted_args);
            
            # Handle varargout in returns
            if ($returns) {
                $returns =~ s/\bvarargout\b/VariableOutputType varargout/g;
            }
            
            return "function $func_name($arguments);";
        }
        # Fallback for simple function signatures
        elsif ($signature =~ /^(\w+)\s*\((.*?)\)\s*$/) {
            my $func_name = $1;
            my $arguments = $2 || '';
            
            # Handle special cases first: varargin
            $arguments =~ s/\bvarargin\b/VariableInputType varargin/g;
            
            # Convert other arguments with InputType, but skip already converted ones
            my @arg_parts = split /\s*,\s*/, $arguments;
            my @converted_args = ();
            
            for my $arg (@arg_parts) {
                $arg =~ s/^\s+|\s+$//g; # trim
                if ($arg eq '') {
                    # Skip empty arguments
                    next;
                } elsif ($arg eq '~') {
                    push @converted_args, 'ignoredArg';
                } elsif ($arg =~ /^VariableInputType/) {
                    # Already converted varargin
                    push @converted_args, $arg;
                } else {
                    # Regular argument
                    push @converted_args, "InputType $arg";
                }
            }
            
            $arguments = join(', ', @converted_args);
            
            return "function $func_name($arguments);";
        }
    }
    
    return "$signature;"; # fallback with semicolon
}

sub generate_cpp_output {
    my ($file_data, $base_name, $package_info) = @_;
    
    my $output = "";
    
    # Parse the docstring
    my $parsed_doc = parse_docstring($file_data->{'docstring'});
    
    # Generate file header comment (from first line of MATLAB docstring)
    my $file_header = generate_file_header_comment($parsed_doc, $base_name, $package_info);
    $output .= $file_header;
    
    # Add blank line to separate file and function comments
    $output .= "\n";
    
    # Add namespace if available
    if ($package_info->{'namespace'}) {
        $output .= "namespace $package_info->{'namespace'} {\n\n";
    }
    
    # Check if this is a class file
    if ($file_data->{'is_class'}) {
        # For classes, add class docstring before class definition
        my $class_comment = generate_function_comment($parsed_doc, $package_info);
        $output .= $class_comment;
        
        # Prepare content without the docstring for convert_class_content
        my $content_without_docstring = prepare_content_without_docstring($file_data);
        
        # Then process the content using original m2cpp.pl logic
        # Pass the content and class signature instead of filename
        my $class_content = convert_class_content($content_without_docstring, $file_data->{'class_signature'});
        $output .= $class_content;
    } else {
        # For functions, use the original function processing
        # Generate function-specific comment (from detailed description)
        my $function_comment = generate_function_comment($parsed_doc, $package_info);
        $output .= $function_comment;
        
        # Convert MATLAB signature to C++
        my $cpp_signature = convert_matlab_signature_to_cpp($file_data->{'function_signature'}, 0);
        $output .= $cpp_signature;
        $output .= "\n";
    }
    
    # Close namespace if opened
    if ($package_info->{'namespace'}) {
        $output .= "\n} // namespace $package_info->{'namespace'}\n";
    }
    
    return $output;
}