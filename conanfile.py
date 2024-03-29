from conans import CMake, ConanFile


class FlagletConan(ConanFile):
    name = "flaglet"
    version = "1.0.0"
    license = "GPL-3"
    url = "https://github.com/astro-informatics/src_flaglet"
    homepage = "https://github.com/astro-informatics/src_flaglet"
    description = "Fast wavelet transform on the ball"
    settings = "os", "arch", "compiler", "build_type"
    topics = ("Physics", "Astrophysics", "Radio Interferometry")
    options = {"fPIC": [True, False]}
    default_options = {"fPIC": True}
    generators = "cmake"
    exports_sources = [
        "src/main/c/*",
        "include/*.h",
        "CMakeLists.txt",
        "cmake/*.cmake",
        "src/test/c/*.c",
        "src/test/c/*.h",
        "src/test/c/CMakeLists.txt",
    ]

    def requirements(self):
        location = "astro-informatics/stable" if self.in_local_cache else "user/testing"
        self.requires(f"s2let/2.2.2@{location}")

    def configure(self):
        if self.settings.compiler == "Visual Studio":
            del self.options.fPIC
        self.options["ssht"].fPIC = self.options.fPIC
        del self.settings.compiler.libcxx

    @property
    def cmake(self):
        if not hasattr(self, "_cmake"):
            self._cmake = CMake(self)
            self._cmake.definitions["tests"] = True
            self._cmake.definitions["conan_deps"] = True
            self._cmake.definitions["fPIC"] = self.options.fPIC
            self._cmake.configure(build_folder="build")
        return self._cmake

    def build(self):
        from pathlib import Path

        path = Path(self.source_folder)
        build = Path(self.source_folder) / "build"
        build.mkdir(exist_ok=True)
        (path / "conanbuildinfo.cmake").rename(path / "build" / "conanbuildinfo.cmake")
        self.cmake.build()
        self.cmake.test()

    def package(self):
        self.cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["flaglet"]
