classdef TestConfigParser < matlab.unittest.TestCase
    % TESTCONFIGPARSER Unit tests for ConfigParser class.
    %
    %   TestConfigParser provides comprehensive test coverage for the
    %   ConfigParser class functionality including configuration parsing,
    %   path management, and error handling.
    %
    % See also:
    %   physics.scheme.ConfigParser

    properties (TestParameter)
        ConfigType = struct( ...
            'Basic', 'env_basic.txt', ...
            'Options', 'env_options.txt', ...
            'Subdirs', 'env_subdirs.txt');
    end

    properties
        TempDir % Temporary directory for test files
        OriginalPath % Original MATLAB path for restoration
    end

    methods (TestMethodSetup)
        function setup(testCase)
            % SETUP Initialize test environment before each test method.

            testCase.TempDir = tempname;
            mkdir(testCase.TempDir);
            testCase.OriginalPath = path;
            testCase.createTestFiles();
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            % TEARDOWN Clean up test environment after each test method.

            restoredefaultpath;
            path(testCase.OriginalPath);

            if isfolder(testCase.TempDir)
                rmdir(testCase.TempDir, 's');
            end
        end
    end

    methods (Test)
        function testConstructor(testCase)
            % TESTCONSTRUCTOR Test ConfigParser constructor.

            parser = physics.scheme.ConfigParser();
            testCase.verifyTrue(isstruct(parser.Config));
            testCase.verifyEmpty(parser.Config.paths);
            testCase.verifyTrue(isstruct(parser.Validators));
        end

        function testBasicParsing(testCase)
            % TESTBASICPARSING Test basic configuration file parsing.

            parser = physics.scheme.ConfigParser();
            configPath = fullfile(testCase.TempDir, 'env_basic.txt');

            parser = parser.parse(configPath);

            testCase.verifyTrue(isstruct(parser.Config));
            testPath = fullfile(testCase.TempDir, 'test_dir');
            testCase.verifyTrue(contains(path, testPath));
            testCase.verifyEqual(length(parser.Config.paths), 1);
        end

        function testOptionParsing(testCase)
            % TESTOPTIONPARSING Test OPTION directive parsing and type
            %   conversion.

            parser = physics.scheme.ConfigParser();
            configPath = fullfile(testCase.TempDir, 'env_options.txt');

            parser = parser.parse(configPath);

            testCase.verifyEqual(parser.Config.testString, 'Hello World');
            testCase.verifyEqual(parser.Config.testNumber, 42);
            testCase.verifyTrue(isnumeric(parser.Config.testNumber));
            testCase.verifyEqual(parser.Config.testLogical, true);
            testCase.verifyTrue(islogical(parser.Config.testLogical));
        end

        function testSubdirParsing(testCase)
            % TESTSUBDIRPARSING Test SUBDIRS directive parsing.

            parser = physics.scheme.ConfigParser();
            configPath = fullfile(testCase.TempDir, 'env_subdirs.txt');

            parser = parser.parse(configPath);

            testDir = fullfile(testCase.TempDir, 'test_dir');
            testCase.verifyTrue(contains(path, testDir));
            testCase.verifyTrue(contains(path, fullfile(testDir, 'subdir1')));
            testCase.verifyGreaterThan(length(parser.Config.paths), 1);
        end

        function testParameterizedConfig(testCase, ConfigType)
            % TESTPARAMETERIZEDCONFIG Test parsing with different configuration
            %   types.

            parser = physics.scheme.ConfigParser();
            configPath = fullfile(testCase.TempDir, ConfigType);

            parser = parser.parse(configPath);

            testCase.verifyTrue(isstruct(parser.Config));

            if strcmp(ConfigType, 'env_options.txt')
                testCase.verifyEqual(parser.Config.testNumber, 42);
            end

            if strcmp(ConfigType, 'env_subdirs.txt')
                subdir1Path = fullfile(testCase.TempDir, 'test_dir', ...
                    'subdir1');
                testCase.verifyTrue(contains(path, subdir1Path));
            end
        end

        function testReset(testCase)
            % TESTRESET Test path reset functionality.

            parser = physics.scheme.ConfigParser();
            configPath = fullfile(testCase.TempDir, 'env_basic.txt');

            parser = parser.parse(configPath);
            testPath = fullfile(testCase.TempDir, 'test_dir');
            testCase.verifyTrue(contains(path, testPath));

            parser = parser.reset();
            testCase.verifyFalse(contains(path, testPath));
            testCase.verifyEmpty(parser.Config.paths);
        end

        function testMissingFile(testCase)
            % TESTMISSINGFILE Test error handling for missing configuration
            %   files.

            parser = physics.scheme.ConfigParser();
            testCase.verifyError(@() parser.parse('nonexistent.txt'), ...
                'MATLAB:validators:mustBeFile');
        end

        function testMultipleParsing(testCase)
            % TESTMULTIPLEPARSING Test parsing multiple configuration files.

            parser = physics.scheme.ConfigParser();

            parser = parser.parse(fullfile(testCase.TempDir, 'env_basic.txt'));
            testPath = fullfile(testCase.TempDir, 'test_dir');
            testCase.verifyTrue(contains(path, testPath));

            parser = parser.parse(fullfile(testCase.TempDir, ...
                'env_options.txt'));

            testCase.verifyGreaterThan(length(parser.Config.paths), 1);
            testCase.verifyEqual(parser.Config.testString, 'Hello World');
            testCase.verifyEqual(parser.Config.testNumber, 42);
        end

        function testComplexOption(testCase)
            % TESTCOMPLEXOPTION Test parsing of complex OPTION values.

            parser = physics.scheme.ConfigParser();

            complexConfigPath = fullfile(testCase.TempDir, 'env_complex.txt');
            fid = fopen(complexConfigPath, 'w');
            fprintf(fid, 'OPTION=testArray=[1,2,3;4,5,6]\n');
            fprintf(fid, 'OPTION=testExpression=sin(pi/4)\n');
            fprintf(fid, 'OPTION=testEquation=5\n');
            fclose(fid);

            parser = parser.parse(complexConfigPath);

            testCase.verifyEqual(parser.Config.testArray, [1, 2, 3; 4, 5, 6]);
            testCase.verifyEqual(parser.Config.testExpression, sin(pi/4), ...
                'AbsTol', 1e-10);
            testCase.verifyEqual(parser.Config.testEquation, 5);
        end

        function testEmptyConfig(testCase)
            % TESTEMPTYCONFIG Test parsing of empty configuration files.

            parser = physics.scheme.ConfigParser();

            emptyConfigPath = fullfile(testCase.TempDir, 'env_empty.txt');
            fid = fopen(emptyConfigPath, 'w');
            fprintf(fid, '# Just a comment\n\n');
            fclose(fid);

            parser = parser.parse(emptyConfigPath);

            testCase.verifyTrue(isstruct(parser.Config));
            % Config always has 'paths' field, so check if only paths exists
            configFields = fieldnames(parser.Config);
            testCase.verifyEqual(length(configFields), 1);
            testCase.verifyEqual(configFields{1}, 'paths');
        end

        function testAddConfig(testCase)
            % TESTADDCONFIG Test adding configuration options with defaults and validators.

            parser = physics.scheme.ConfigParser();
            
            % Test with default value only
            parser = parser.addConfig('testDefault', default=42);
            testCase.verifyEqual(parser.Config.testDefault, 42);
            
            % Test with validator only
            validator = @(x) x > 0;
            parser = parser.addConfig('testValidator', validator=validator);
            testCase.verifyEqual(parser.Validators.testValidator, validator);
            
            % Test with both default and validator
            parser = parser.addConfig('testBoth', default=10, validator=@(x) x > 5);
            testCase.verifyEqual(parser.Config.testBoth, 10);
            testCase.verifyTrue(isa(parser.Validators.testBoth, 'function_handle'));
        end

        function testGetConfig(testCase)
            % TESTGETCONFIG Test getting configuration option values.

            parser = physics.scheme.ConfigParser();
            parser = parser.addConfig('testOption', default='testValue');
            
            % Test successful retrieval
            value = parser.getConfig('testOption');
            testCase.verifyEqual(value, 'testValue');
            
            % Test error for non-existent option
            testCase.verifyError(@() parser.getConfig('nonExistent'), ...
                'physics:scheme:ConfigParser:ConfigNotFound');
        end

        function testSetConfig(testCase)
            % TESTSETCONFIG Test setting configuration option values with validation.

            parser = physics.scheme.ConfigParser();
            
            % Test setting without validator
            parser = parser.addConfig('testNoValidator', default=0);
            parser = parser.setConfig('testNoValidator', 42);
            testCase.verifyEqual(parser.Config.testNoValidator, 42);
            
            % Test setting with validator (valid value)
            parser = parser.addConfig('testWithValidator', default=5, validator=@(x) x > 0);
            parser = parser.setConfig('testWithValidator', 10);
            testCase.verifyEqual(parser.Config.testWithValidator, 10);
            
            % Test validation failure
            testCase.verifyError(@() parser.setConfig('testWithValidator', -5), ...
                'physics:scheme:ConfigParser:InvalidConfig');
        end
    end

    methods (Access = private)
        function createTestFiles(testCase)
            % CREATETESTFILES Create test configuration files and directories.

            testDir = fullfile(testCase.TempDir, 'test_dir');
            mkdir(testDir);
            mkdir(fullfile(testDir, 'subdir1'));
            mkdir(fullfile(testDir, 'subdir2'));

            basicConfig = fopen(fullfile(testCase.TempDir, ...
                'env_basic.txt'), 'w');
            fprintf(basicConfig, 'MATLABPATH=%s\n', testDir);
            fclose(basicConfig);

            optionsConfig = fopen(fullfile(testCase.TempDir, ...
                'env_options.txt'), 'w');
            fprintf(optionsConfig, 'MATLABPATH=%s\n', testDir);
            fprintf(optionsConfig, 'OPTION=testString=Hello World\n');
            fprintf(optionsConfig, 'OPTION=testNumber=42\n');
            fprintf(optionsConfig, 'OPTION=testLogical=true\n');
            fclose(optionsConfig);

            subdirsConfig = fopen(fullfile(testCase.TempDir, ...
                'env_subdirs.txt'), 'w');
            fprintf(subdirsConfig, 'SUBDIRS=%s\n', testDir);
            fclose(subdirsConfig);
        end
    end
end