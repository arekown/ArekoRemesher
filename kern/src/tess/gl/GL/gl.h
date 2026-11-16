/* Types and constants the polygon tessellator sources need; nothing else of
 * a graphics library is used by the core. */
#ifndef ARM_TESS_GL_H
#define ARM_TESS_GL_H

typedef unsigned int GLenum;
typedef unsigned char GLboolean;
typedef unsigned char GLubyte;
typedef int GLint;
typedef int GLsizei;
typedef float GLfloat;
typedef double GLdouble;
typedef void GLvoid;

#define GL_FALSE 0
#define GL_TRUE 1
#define GL_LINE_LOOP 0x0002
#define GL_TRIANGLES 0x0004
#define GL_TRIANGLE_STRIP 0x0005
#define GL_TRIANGLE_FAN 0x0006

#endif
