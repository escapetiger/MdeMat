classdef TestSchemeHierarchy < matlab.unittest.TestCase
    % TESTSCHEMEHIERARCHY Test the redesigned Scheme hierarchy.
    %
    %   TestSchemeHierarchy verifies that the new Scheme hierarchy properly
    %   inherits from ConfigParser and provides the expected configuration
    %   management capabilities.
    %
    % See also:
    %   physics.scheme.Scheme, physics.scheme.ConfigParser

    methods (Test)
        function testSchemeInheritance(testCase)
            % TESTSCHEMEINHERITANCE Test that Scheme inherits from ConfigParser.
            
            scheme = physics.scheme.Scheme();
            testCase.verifyTrue(isa(scheme, 'physics.scheme.ConfigParser'));
            testCase.verifyTrue(isa(scheme, 'physics.scheme.Scheme'));
        end

        function testSchemeDefaultConfiguration(testCase)
            % TESTSCHEMEDEFAULTCONFIGURATION Test default configuration setup.
            
            scheme = physics.scheme.Scheme();
            
            % Verify default verbose setting
            testCase.verifyEqual(scheme.getConfig('verbose'), 0);
            
            % Verify Config structure exists
            testCase.verifyTrue(isstruct(scheme.Config));
            testCase.verifyTrue(isfield(scheme.Config, 'verbose'));
        end

        function testSchemeStructConfiguration(testCase)
            % TESTSCHEMESTRUCTCONFIGURATION Test configuration from struct.
            
            config = struct('verbose', 2, 'customParam', 'testValue');
            scheme = physics.scheme.Scheme(config=config);
            
            % Verify configured values
            testCase.verifyEqual(scheme.getConfig('verbose'), 2);
            testCase.verifyEqual(scheme.getConfig('customParam'), 'testValue');
        end

        function testOdeSchemeInheritance(testCase)
            % TESTODESCHEMEINHERTIANCE Test OdeScheme inheritance chain.
            
            scheme = physics.scheme.OdeScheme();
            testCase.verifyTrue(isa(scheme, 'physics.scheme.ConfigParser'));
            testCase.verifyTrue(isa(scheme, 'physics.scheme.Scheme'));
            testCase.verifyTrue(isa(scheme, 'physics.scheme.OdeScheme'));
        end

        function testOdeSchemeDefaults(testCase)
            % TESTODESCHEMEDEFAULTS Test OdeScheme default configuration.
            
            scheme = physics.scheme.OdeScheme();
            
            % Verify parent defaults
            testCase.verifyEqual(scheme.getConfig('verbose'), 0);
            
            % Verify ODE-specific defaults
            testCase.verifyEqual(scheme.getConfig('tOdeIntName'), 'be');
            testCase.verifyEqual(scheme.getConfig('tFinal'), 1.0);
            testCase.verifyEmpty(scheme.getConfig('tSplitName'));
        end

        function testOdeSchemeConfiguration(testCase)
            % TESTODESCHEMECONFIGURATTION Test OdeScheme configuration.
            
            config = struct('tFinal', 10.0, 'tOdeIntName', 'sdirk3');
            scheme = physics.scheme.OdeScheme(config=config);
            
            testCase.verifyEqual(scheme.getConfig('tFinal'), 10.0);
            testCase.verifyEqual(scheme.Config.tOdeIntName, {'sdirk3'});
        end

        function testOdeSchemeValidation(testCase)
            % TESTODESCHEMEVALIDATION Test OdeScheme validation.
            
            % Test invalid tFinal
            config = struct('tFinal', -1.0);
            testCase.verifyError(@() physics.scheme.OdeScheme(config=config), ...
                'physics:scheme:ConfigParser:InvalidConfig');
            
            % Test invalid tOdeIntName
            config = struct('tOdeIntName', 'invalid_method');
            testCase.verifyError(@() physics.scheme.OdeScheme(config=config), ...
                'physics:scheme:ConfigParser:InvalidConfig');
        end

        function testConfigFileCreation(testCase)
            % TESTCONFIGFILECREATION Test configuration file parsing capability.
            
            % Create temporary config file
            tempFile = tempname;
            fid = fopen(tempFile, 'w');
            fprintf(fid, 'OPTION=verbose=1\n');
            fprintf(fid, 'OPTION=tFinal=5.0\n');
            fprintf(fid, 'OPTION=tOdeIntName=bdf2\n');
            fclose(fid);
            
            % Test file-based configuration
            scheme = physics.scheme.OdeScheme(file=tempFile);
            
            testCase.verifyEqual(scheme.getConfig('verbose'), 1);
            testCase.verifyEqual(scheme.getConfig('tFinal'), 5.0);
            testCase.verifyEqual(scheme.Config.tOdeIntName, {'bdf2'});
            
            % Cleanup
            delete(tempFile);
        end

        function testMixedConfiguration(testCase)
            % TESTMIXEDCONFIGURATION Test file + struct configuration.
            
            % Create temporary config file
            tempFile = tempname;
            fid = fopen(tempFile, 'w');
            fprintf(fid, 'OPTION=verbose=1\n');
            fprintf(fid, 'OPTION=tFinal=5.0\n');
            fclose(fid);
            
            % Override with struct
            config = struct('verbose', 2, 'tOdeIntName', 'sdirk2');
            scheme = physics.scheme.OdeScheme(file=tempFile, config=config);
            
            % Struct should override file values
            testCase.verifyEqual(scheme.getConfig('verbose'), 2);
            testCase.verifyEqual(scheme.getConfig('tFinal'), 5.0);
            testCase.verifyEqual(scheme.Config.tOdeIntName, {'sdirk2'});
            
            % Cleanup
            delete(tempFile);
        end
    end
end