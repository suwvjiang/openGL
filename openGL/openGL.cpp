// openGL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include <Windows.h>
//#include <GL/glew.h>
#include <GL/glut.h>

//#pragma comment(lib, "glew32.lib")

void init(void)
{
	glClearColor(1.0f, 1.0f, 0.0f, 0.0f);//设定背景颜色

	glMatrixMode(GL_PROJECTION);//设定投影模式
	gluOrtho2D(0.0, 200.0, 0, 150.0);
}

void lineSegment(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.0f, 0.4f, 0.2f);//设定显示对象的颜色
	glBegin(GL_LINES);
	glVertex2i(180, 15);
	glVertex2i(10, 145);
	glEnd();

	glFlush();
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);//初始化

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);//设定窗口的缓存和颜色模型

	glutInitWindowPosition(100, 100);//设定窗口位置

	glutInitWindowSize(800, 600);//设定窗口大小

	glutCreateWindow("Tutorial 01");//设置窗口标题

	init();
	glutDisplayFunc(lineSegment);//指定显示内容
	glutMainLoop();
	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
