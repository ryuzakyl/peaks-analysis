/*
    Copyright (C) CENATAV, DATYS - All Rights Reserved
    Unauthorized copying of this file, via any medium is strictly prohibited
    Proprietary and confidential
    Written by Victor M. Mendiola Lau <vmendiola@cenatav.co.cu>, March 2017
*/

#ifndef CONFIG_HPP_INCLUDED
#define CONFIG_HPP_INCLUDED

/*
	cross platform exportation logic definitions
*/

#ifndef EXTERN_C
#  ifdef __cplusplus
#    define EXTERN_C extern "C"
#  else
#    define EXTERN_C
#  endif
#endif

#if defined WIN32 || defined _WIN32
#  define CDECL __cdecl
#  define STDCALL __stdcall
#else
#  define CDECL
#  define STDCALL
#endif

#if defined WIN32 || defined _WIN32
#  define EXPORTS __declspec(dllexport)
#else
#  define EXPORTS
#endif

#ifndef API_FUNC
#  define API_FUNC(ret_type) EXTERN_C EXPORTS ret_type CDECL
#endif

#endif // CONFIG_HPP_INCLUDED
