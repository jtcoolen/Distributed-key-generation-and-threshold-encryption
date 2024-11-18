with import <nixpkgs> { };
llvmPackages.libcxxStdenv.mkDerivation {
  name = "env";
  nativeBuildInputs = [ clang-tools clang ];
  buildInputs = [
    clang
    gmp
    pkg-config
    m4
    cmake
    libev
    libiconv
    autoconf
    gnumake
    cmake
    openssl
    bison
    gnumake
    pkg-config
    gcc
    gmp
    gmpxx
    libcxx
    cargo
    rustc
    ninja
  ] ++ lib.optionals libcxxStdenv.isDarwin
    (with darwin.apple_sdk.frameworks; [ Cocoa CoreServices Security ]);

  NIX_LDFLAGS =
    # Work around https://github.com/NixOS/nixpkgs/issues/166205.
    lib.optionalString (stdenv.cc.isClang && stdenv.cc.libcxx != null)
    " -l${stdenv.cc.libcxx.cxxabi.libName} -lc -lm -lgmp /opt/homebrew/lib/libgmp.dylib /opt/homebrew/lib/libgmpxx.dylib /opt/homebrew/include/gmp.h";

  shellHook = ''
    export NIX_LDFLAGS="-F${pkgs.darwin.apple_sdk.frameworks.CoreFoundation}/Library/Frameworks -framework CoreFoundation $NIX_LDFLAGS";
  '';
}
