with import <nixpkgs> { };
llvmPackages.libcxxStdenv.mkDerivation {
  name = "env";
  nativeBuildInputs = [ clang clang-tools ];
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
    gmp
    gmpxx
    libcxx
    cargo
    rustc
    ninja
    valgrind
  ];
 }
