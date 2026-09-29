# Tensorlab (required for CPD / tensor steps)

Tensorlab files are included in this folder for the thesis CPD / tensor
experiments. Tensorlab is distributed by its authors via
[tensorlab.net](http://www.tensorlab.net/versions.html); see
`license_tensorlab.txt` for the bundled license text.

## MATLAB path setup

From the project root, add Tensorlab recursively:

```matlab
addpath(genpath(fullfile(pwd, 'tensorlab')));
savepath;  % optional
```

The SSA scripts in this project **do not** require Tensorlab; only the later thesis tensor / CPD experiments do.

## Verify

After adding the path, run Tensorlab’s own demos or:

```matlab
which cpd_nls -all
```

If that resolves to a file under your Tensorlab folder, the path is set correctly.
