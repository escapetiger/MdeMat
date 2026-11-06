function root = toolboxRoot()
%TOOLBOXROOT Locate distributable resources independently of caller pwd.
root = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));
end
