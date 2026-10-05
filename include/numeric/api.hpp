#ifndef ARIADNE_NUMERIC_API_HPP
#define ARIADNE_NUMERIC_API_HPP

#if defined(_WIN32)
  #if defined(ARIADNE_NUMERIC_STATIC) || defined(ARIADNE_NUMERIC_EMBEDDED)
    #define ARIADNE_NUMERIC_API
  #elif defined(ARIADNE_NUMERIC_BUILD)
    #define ARIADNE_NUMERIC_API __declspec(dllexport)
  #else
    #define ARIADNE_NUMERIC_API __declspec(dllimport)
  #endif
#else
  #define ARIADNE_NUMERIC_API
#endif

#endif
