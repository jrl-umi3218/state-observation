{
  description = "Nix flake for state-obervation";

  inputs.mc-rtc-nix.url = "github:mc-rtc/nixpkgs";

  outputs =
    inputs:
    inputs.mc-rtc-nix.lib.mkFlakoboros inputs (
      { lib, ... }:
      {
        overrideAttrs.state-observation =
          { drv-prev, pkgs-final, ... }:
          {
            src = lib.cleanSource ./.;
            nativeBuildInputs = [
              pkgs-final.doxygen
              pkgs-final.graphviz
              pkgs-final.gbenchmark
              pkgs-final.jrl-cmakemodulesv2
            ]
            ++ drv-prev.nativeBuildInputs;
            cmakeFlags = [ (lib.cmakeBool "INSTALL_DOCUMENTATION" true) ] ++ drv-prev.cmakeFlags;
            doCheck = true;
          };
      }
    );
}
