#ifndef DEBUG_H
#define DEBUG_H

#define DEBUG_MODE // Comment/uncomment this line to toggle debug mode
// #define TRACE_DEBUG_MODE // Comment/uncomment for very verbose debug messages

#if defined DEBUG_MODE || defined TRACE_DEBUG_MODE
#define DEBUG(x) x
#else 
#define DEBUG(x)
#endif

#ifdef TRACE_DEBUG_MODE
#define TRACE(x) x
#else 
#define TRACE(x)
#endif

#endif
