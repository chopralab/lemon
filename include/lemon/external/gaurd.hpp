#ifndef LEMON_EXTERNAL_GUARD_HPP

#ifdef _MSVC_LANG

#define LEMON_EXTERNAL_FILE_PUSH \
	__pragma(warning(push)) \
	__pragma(warning(disable : 4191)) \
	__pragma(warning(disable : 4242)) \
	__pragma(warning(disable : 4365)) \
	__pragma(warning(disable : 4464)) \
	__pragma(warning(disable : 4466)) \
	__pragma(warning(disable : 4505)) \
	__pragma(warning(disable : 4643)) \
	__pragma(warning(disable : 4868)) \
    __pragma(warning(disable : 5039)) \
	__pragma(warning(disable : 26437)) \
	__pragma(warning(disable : 26444)) \
	__pragma(warning(disable : 26451)) \
	__pragma(warning(disable : 26495)) \
	__pragma(warning(disable : 26812)) \
    __pragma(warning(disable : 28182))

#define LEMON_EXTERNAL_FILE_POP	__pragma(warning(pop))

#else

#define LEMON_EXTERNAL_FILE_PUSH
#define LEMON_EXTERNAL_FILE_POP

#endif // _MSVC_LANG

#endif // !LEMON_EXTERNAL_GUARD_HPP
