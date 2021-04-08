#ifndef FLAGLET_VERSION_H
#define FLAGLET_VERSION_H
inline const char *flaglet_version_string() { return "@PROJECT_VERSION@"; }
inline const char *flaglet_info() {
  return "package:\n"
         "  name: FLAGLET\n"
         "  description: Fast wavelet transform on the ball\n"
         "  authors:\n"
         "      - Boris Leistedt\n"
         "      - Jason McEwen\n"
         "  license: GPL-3\n"
         "  url: https://github.com/astro-informatics/flaglet\n"
         "  version: @PROJECT_VERSION@\n";
};
// clang-format off
inline int flaglet_version_major() { return @PROJECT_VERSION_MAJOR@; }
inline int flaglet_version_minor() { return @PROJECT_VERSION_MINOR@; }
inline int flaglet_version_patch() { return @PROJECT_VERSION_PATCH@; }
#define FLAGLET_VERSION_MAJOR  @PROJECT_VERSION_MAJOR@
#define FLAGLET_VERSION_MINOR  @PROJECT_VERSION_MINOR@
#define FLAGLET_VERSION_PATCH  @PROJECT_VERSION_PATCH@
// clang-format on
#endif
