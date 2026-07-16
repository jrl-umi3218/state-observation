{
  description = "Describes interfaces for state observers, and implements some observers (including linear and extended Kalman filters)";

  inputs.mc-rtc-nix.url = "github:mc-rtc/nixpkgs";

  outputs =
    inputs:
    inputs.mc-rtc-nix.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        overrideAttrs.state-observation =
          { drv-prev, pkgs-final, ... }:
          {
            outputs = [
              "out"
              "doc"
            ];
            src = lib.cleanSource ./.;
            nativeBuildInputs = [
              pkgs-final.doxygen
              pkgs-final.graphviz
              pkgs-final.gbenchmark
              pkgs-final.jrl-cmakemodulesv2
            ]
            ++ drv-prev.nativeBuildInputs;
            cmakeFlags = drv-prev.cmakeFlags ++ [ (lib.cmakeBool "INSTALL_DOCUMENTATION" true) ];
            doCheck = true;
          };
      }
    );
}
