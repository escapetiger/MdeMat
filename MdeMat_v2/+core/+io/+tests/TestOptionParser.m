classdef TestOptionParser < matlab.unittest.TestCase
    % TESTOPTIONPARSER Unit tests for OptionParser class.
    %
    %   TestOptionParser provides comprehensive test coverage for the
    %   OptionParser class functionality including configuration parsing,
    %   path management, and error handling.
    %
    % Examples:
    %   % Run all tests
    %   results = run(core.io.tests.TestOptionParser);
    %
    % See Also:
    %   core.io.OptionParser

    properties (TestParameter)
        configType = struct( ...
            'Basic', 'env_basic.txt', ...
            'Options', 'env_options.txt', ...
            'Subdirs', 'env_subdirs.txt');
    end

    properties
        tempDir % Temporary directory for test files
        originalPath % Original MATLAB path for restoration
    end

    methods (TestMethodSetup)
        function setup(testCase)
            % SETUP Initialize test environment before each test method.

            testCase.tempDir = tempname;
            mkdir(testCase.tempDir);
            testCase.originalPath = path;
            testCase.createTestFiles();
        end
    end

    methods (TestMethodTeardown)
        function teardown(testCase)
            % TEARDOWN Clean up test environment after each test method.

            restoredefaultpath;
            path(testCase.originalPath);

            if isfolder(testCase.tempDir)
                rmdir(testCase.tempDir, 's');
            end
        end
    end

    methods (Test)
        function testConstructor(testCase)
            % TESTCONSTRUCTOR Test OptionParser constructor.

            parser = core.io.OptionParser();
            testCase.verifyEmpty(parser.paths);
        end

        function testBasicParsing(testCase)
            % TESTBASICPARSING Test basic configuration file parsing.

            parser = core.io.OptionParser();
            configPath = fullfile(testCase.tempDir, 'env_basic.txt');

            options = parser.parse(configPath);

            testCase.verifyTrue(isstruct(options));
            testCase.verifyTrue(contains(path, fullfile(testCase.tempDir, 'test_dir')));
            testCase.verifyEqual(length(parser.paths), 1);
        end

        function testOptionParsing(testCase)
            % TESTOPTIONPARSING Test OPTION directive parsing and type conversion.

            parser = core.io.OptionParser();
            configPath = fullfile(testCase.tempDir, 'env_options.txt');

            options = parser.parse(configPath);

            testCase.verifyEqual(options.testString, 'Hello World');
            testCase.verifyEqual(options.testNumber, 42);
            testCase.verifyTrue(isnumeric(options.testNumber));
            testCase.verifyEqual(options.testLogical, true);
            testCase.verifyTrue(islogical(options.testLogical));
        end

        function testSubdirParsing(testCase)
            % TESTSUBDIRPARSING Test SUBDIRS directive parsing.

            parser = core.io.OptionParser();
            configPath = fullfile(testCase.tempDir, 'env_subdirs.txt');

            parser.parse(configPath);

            testCase.verifyTrue(contains(path, fullfile(testCase.tempDir, 'test_dir')));
            testCase.verifyTrue(contains(path, fullfile(testCase.tempDir, 'test_dir', 'subdir1')));
            testCase.verifyGreaterThan(length(parser.paths), 1);
        end

        function testParameterizedConfig(testCase, configType)
            % TESTPARAMETERIZEDCONFIG Test parsing with different configuration types.

            parser = core.io.OptionParser();
            configPath = fullfile(testCase.tempDir, configType);

            options = parser.parse(configPath);

            testCase.verifyTrue(isstruct(options));

            if strcmp(configType, 'env_options.txt')
                testCase.verifyEqual(options.testNumber, 42);
            end

            if strcmp(configType, 'env_subdirs.txt')
                testCase.verifyTrue(contains(path, fullfile(testCase.tempDir, 'test_dir', 'subdir1')));
            end
        end

        function testReset(testCase)
            % TESTRESET Test path reset functionality.

            parser = core.io.OptionParser();
            configPath = fullfile(testCase.tempDir, 'env_basic.txt');

            parser.parse(configPath);
            testCase.verifyTrue(contains(path, fullfile(testCase.tempDir, 'test_dir')));

            parser.reset();
            testCase.verifyFalse(contains(path, fullfile(testCase.tempDir, 'test_dir')));
            testCase.verifyEmpty(parser.paths);
        end

        function testMissingFile(testCase)
            % TESTMISSINGFILE Test error handling for missing configuration files.

            parser = core.io.OptionParser();
            testCase.verifyError(@() parser.parse('nonexistent.txt'), ...
                'core:io:OptionParser:FileNotFound');
        end

        function testMultipleParsing(testCase)
            % TESTMULTIPLEPARSING Test parsing multiple configuration files.

            parser = core.io.OptionParser();

            parser.parse(fullfile(testCase.tempDir, 'env_basic.txt'));
            testCase.verifyTrue(contains(path, fullfile(testCase.tempDir, 'test_dir')));

            options = parser.parse(fullfile(testCase.tempDir, 'env_options.txt'));

            testCase.verifyGreaterThan(length(parser.paths), 1);
            testCase.verifyEqual(options.testString, 'Hello World');
            testCase.verifyEqual(options.testNumber, 42);
        end

        function testComplexOption(testCase)
            % TESTCOMPLEXOPTION Test parsing of complex OPTION values.

            parser = core.io.OptionParser();

            complexConfigPath = fullfile(testCase.tempDir, 'env_complex.txt');
            fid = fopen(complexConfigPath, 'w');
            fprintf(fid, 'OPTION=testArray=[1,2,3;4,5,6]\n');
            fprintf(fid, 'OPTION=testExpression=sin(pi/4)\n');
            fprintf(fid, 'OPTION=testEquation=5\n');
            fclose(fid);

            options = parser.parse(complexConfigPath);

            testCase.verifyEqual(options.testArray, [1, 2, 3; 4, 5, 6]);
            testCase.verifyEqual(options.testExpression, sin(pi/4), 'AbsTol', 1e-10);
            testCase.verifyEqual(options.testEquation, 5);
        end

        function testInvalidConfig(testCase)
            % TESTINVALIDCONFIG Test warning handling for invalid configuration lines.

            parser = core.io.OptionParser();

            InvalidConfigPath = fullfile(testCase.tempDir, 'env_invalid.txt');
            fid = fopen(InvalidConfigPath, 'w');
            fprintf(fid, 'INVALID_LINE_NO_EQUALS\n');
            fprintf(fid, '# This is a comment\n');
            fprintf(fid, 'OPTION=validOption=123\n');
            fclose(fid);

            testCase.verifyWarning(@() parser.parse(InvalidConfigPath), ...
                'core:io:OptionParser:InvalidConfig');
        end

        function testEmptyConfig(testCase)
            % TESTEMPTYCONFIG Test parsing of empty configuration files.

            parser = core.io.OptionParser();

            emptyConfigPath = fullfile(testCase.tempDir, 'env_empty.txt');
            fid = fopen(emptyConfigPath, 'w');
            fprintf(fid, '# Just a comment\n\n');
            fclose(fid);

            options = parser.parse(emptyConfigPath);

            testCase.verifyTrue(isstruct(options));
            testCase.verifyTrue(isempty(fieldnames(options)));
        end
    end

    methods (Access = private)
        function createTestFiles(testCase)
            % CREATETESTFILES Create test configuration files and directories.

            testDir = fullfile(testCase.tempDir, 'test_dir');
            mkdir(testDir);
            mkdir(fullfile(testDir, 'subdir1'));
            mkdir(fullfile(testDir, 'subdir2'));

            basicConfig = fopen(fullfile(testCase.tempDir, 'env_basic.txt'), 'w');
            fprintf(basicConfig, 'MATLABPATH=%s\n', testDir);
            fclose(basicConfig);

            optionsConfig = fopen(fullfile(testCase.tempDir, 'env_options.txt'), 'w');
            fprintf(optionsConfig, 'MATLABPATH=%s\n', testDir);
            fprintf(optionsConfig, 'OPTION=testString=Hello World\n');
            fprintf(optionsConfig, 'OPTION=testNumber=42\n');
            fprintf(optionsConfig, 'OPTION=testLogical=true\n');
            fclose(optionsConfig);

            subdirsConfig = fopen(fullfile(testCase.tempDir, 'env_subdirs.txt'), 'w');
            fprintf(subdirsConfig, 'SUBDIRS=%s\n', testDir);
            fclose(subdirsConfig);
        end
    end
end