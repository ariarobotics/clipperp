/*
Using the debug macro in the code: 
Surround the debugging code with #ifdef and #endif directives based
on the state of the DEBUG macro. This will include or exclude the 
debugging code depending on whether the macro is defined.

Compile with or without the -D flag: 
To enable or disable debugging, you can compile your code with or 
without the -D flag, which defines or undefines the DEBUG macro, respectively. 

*/

// debug.h
#ifndef DEBUG_H
#define DEBUG_H

// Uncomment the line below to enable debugging
// #define DEBUG
#define DEBUG_TIMING

#endif // DEBUG_H

