function description = describe(options)
%DESCRIBE Describe this native MATLAB process without personal/license data.
%   installedProducts=true explicitly adds ver's installation inventory.
%   Installed products do not establish runtime dependencies or available
%   licenses. No hostname, username, license number, timestamp, seed or Git
%   version is inferred. Supply execution provenance explicitly to RunRecord.
arguments
    options.installedProducts (1, 1) logical = false
end
identity = mdemat.tools.environment.versionInfo();
products = struct('Name', {}, 'Version', {}, 'Release', {}, 'Date', {});
if options.installedProducts, products = ver; end
description = struct('Type', 'mdemat.environment', 'SchemaVersion', 1, ...
    'MatlabVersion', version, 'MatlabRelease', version('-release'), 'Computer', computer, ...
    'ToolboxVersion', identity.ToolboxVersion, ...
    'InstalledProductsIncluded', options.installedProducts, 'InstalledProducts', products, ...
    'ProductMeaning', 'Installation inventory only; no dependency or license-availability inference.');
end
