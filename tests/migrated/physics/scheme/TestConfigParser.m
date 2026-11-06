classdef TestConfigParser < matlab.unittest.TestCase
    % Pure readConfig preserves file parsing while separating resource paths.
    properties
        Folder
    end
    methods (TestMethodSetup)
        function temporaryDirectory(test)
            fixture = test.applyFixture(matlab.unittest.fixtures.TemporaryFolderFixture);
            test.Folder = fixture.Folder;
        end
    end
    methods (Test)
        function scalarAndRectangularValuesRetainTheirTypes(test)
            file = test.write({'OPTION=text=''Hello World''', 'OPTION=number=42', ...
                'OPTION=enabled=true', 'OPTION=matrix=[1,2,3;4,5,6]', ...
                'OPTION=cells={''all'', ''right''; ''left'', ''''}', ...
                'OPTION=label=''a=b#c'' # real comment', 'OPTION=empty=[]'});
            original = path;
            [values, info] = mdemat.tools.environment.readConfig(file);
            test.verifyEqual(values.text, 'Hello World');
            test.verifyEqual(values.number, 42);
            test.verifyTrue(values.enabled);
            test.verifyClass(values.enabled, 'logical');
            test.verifyEqual(values.matrix, [1,2,3;4,5,6]);
            test.verifyEqual(values.cells, {'all','right';'left',''});
            test.verifyEqual(values.label, 'a=b#c');
            test.verifyEmpty(values.empty);
            test.verifyEqual(info.Lines.number, 2);
            test.verifyEqual(path, original);
        end

        function arithmeticMatchesIndependentMatlabOperatorValues(test)
            file = test.write({'OPTION=a=2^3^2', 'OPTION=b=-2^2', ...
                'OPTION=c=2^-2', 'OPTION=d=2^-3^2', 'OPTION=e=-2^-3^2', ...
                'OPTION=f=2^+3^2', 'OPTION=g=sin(pi/4)', ...
                'OPTION=h=2 / 3', 'OPTION=bracketed=[(1 + 2), 4]'});
            values = mdemat.tools.environment.readConfig(file);
            test.verifyEqual([values.a, values.b, values.c, values.d, values.e, values.f], ...
                [64, -4, 0.25, 0.015625, -0.015625, 64]);
            test.verifyEqual(values.g, 1/sqrt(2), AbsTol=8*eps);
            test.verifyEqual(values.h, 2/3);
            test.verifyEqual(values.bracketed, [3, 4]);
        end

        function pathsAreDescriptionsSeparateFromOptionNames(test)
            mkdir(fullfile(test.Folder, 'shared', 'nested'));
            file = test.write({'OPTION=paths=17', 'MATLABPATH=shared', 'SUBDIRS=shared'});
            original = path;
            [values, info] = mdemat.tools.environment.readConfig(file);
            test.verifyEqual(values.paths, 17);
            test.verifyTrue(ismember(fullfile(test.Folder, 'shared'), info.Paths));
            test.verifyTrue(ismember(fullfile(test.Folder, 'shared', 'nested'), info.Paths));
            test.verifyEqual(numel(info.Paths), 2);
            test.verifyEqual(path, original);
        end

        function relativeResourcesDoNotBindToCallerDirectory(test)
            mkdir(fullfile(test.Folder, 'shared'));
            caller = fullfile(test.Folder, 'caller');
            mkdir(fullfile(caller, 'shared'));
            file = test.write({'MATLABPATH=shared'});
            original = pwd;
            guard = onCleanup(@() cd(original)); %#ok<NASGU>
            cd(caller);
            [~, info] = mdemat.tools.environment.readConfig(file);
            test.verifyEqual(info.Paths, {fullfile(test.Folder, 'shared')});
            test.verifyEqual(pwd, caller);
        end

        function unicodeBomAndEscapedQuotesAreText(test)
            file = test.write({'OPTION=label="中文 路径"', ...
                'OPTION=quote=''can''''t # stop''', 'OPTION=literal="a=b#c"'});
            stream = fopen(file, 'r', 'n', 'UTF-8');
            text = fscanf(stream, '%c');
            fclose(stream);
            stream = fopen(file, 'w', 'n', 'UTF-8');
            fprintf(stream, '%s%s', char(65279), text);
            fclose(stream);
            values = mdemat.tools.environment.readConfig(file);
            test.verifyEqual(values.label, "中文 路径");
            test.verifyEqual(values.quote, 'can''t # stop');
            test.verifyEqual(values.literal, "a=b#c");
        end

        function malformedAndExecutableExpressionsAreRejected(test)
            expressions = {'feval(''sin'', pi)', '[1,,2]', '[1,]', ...
                '[;1]', '[(1+2}]', '1 e2', 'system(''anything'')'};
            for index = 1:numel(expressions)
                file = test.write({['OPTION=value=', expressions{index}]});
                test.verifyError(@() mdemat.tools.environment.readConfig(file), ...
                    'mdemat:tools:environment:readConfig:UnsupportedExpression');
            end
        end

        function duplicateUnknownAndMalformedDirectivesHaveContext(test)
            file = test.write({'OPTION=value=1', 'OPTION=value=2'});
            test.verifyError(@() mdemat.tools.environment.readConfig(file), ...
                'mdemat:tools:environment:readConfig:DuplicateOption');
            file = test.write({'# first line', 'BAD=value'});
            try
                mdemat.tools.environment.readConfig(file);
                test.assertFail('Expected an unknown directive to fail.');
            catch exception
                test.verifyEqual(exception.identifier, ...
                    'mdemat:tools:environment:readConfig:UnknownSetting');
                test.verifySubstring(exception.message, [file, ':2:']);
                test.verifyNotEmpty(exception.cause);
            end
            file = test.write({'OPTION=bad name=1'});
            test.verifyError(@() mdemat.tools.environment.readConfig(file), ...
                'mdemat:tools:environment:readConfig:InvalidOption');
        end

        function emptyIndependentAndMissingFiles(test)
            file = test.write({'# comment', ''});
            [values, info] = mdemat.tools.environment.readConfig(file);
            test.verifyEmpty(fieldnames(values));
            test.verifyEmpty(info.Paths);
            second = test.write({'OPTION=value=42'});
            parsed = mdemat.tools.environment.readConfig(second);
            test.verifyEqual(parsed.value, 42);
            test.verifyEmpty(fieldnames(mdemat.tools.environment.readConfig(file)));
            test.verifyError(@() mdemat.tools.environment.readConfig( ...
                fullfile(test.Folder, 'missing.txt')), 'MATLAB:validators:mustBeFile');
        end

        function normalEofIncludingNoFinalNewlineIsNotAnIoFailure(test)
            payloads = {'', 'OPTION=value=1', sprintf('OPTION=value=2\n')};
            for index = 1:numel(payloads)
                file = [tempname(test.Folder), '.txt'];
                stream = fopen(file, 'w', 'n', 'UTF-8');
                fprintf(stream, '%s', payloads{index});
                fclose(stream);
                values = mdemat.tools.environment.readConfig(file);
                if index == 1
                    test.verifyEmpty(fieldnames(values));
                else
                    test.verifyEqual(values.value, index-1);
                end
            end
        end
    end
    methods (Access = private)
        function file = write(test, lines)
            file = [tempname(test.Folder), '.txt'];
            stream = fopen(file, 'w', 'n', 'UTF-8');
            test.assertGreaterThanOrEqual(stream, 0);
            guard = onCleanup(@() fclose(stream)); %#ok<NASGU>
            for index = 1:numel(lines), fprintf(stream, '%s\n', lines{index}); end
        end
    end
end
