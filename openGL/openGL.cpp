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
	glClearColor(1.0f, 1.0f, 0.0f, 0.0f);//�趨������ɫ

	glMatrixMode(GL_PROJECTION);//�趨ͶӰģʽ
	gluOrtho2D(0.0, 200.0, 0, 150.0);
}

void lineSegment(void)
{
	glClear(GL_COLOR_BUFFER_BIT);

	glColor3f(0.0f, 0.4f, 0.2f);//�趨��ʾ�������ɫ
	glBegin(GL_LINES);
	glVertex2i(180, 15);
	glVertex2i(10, 145);
	glEnd();

	glFlush();
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);//��ʼ��

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);//�趨���ڵĻ������ɫģ��

	glutInitWindowPosition(100, 100);//�趨����λ��

	glutInitWindowSize(800, 600);//�趨���ڴ�С

	glutCreateWindow("Tutorial 01");//���ô��ڱ���

	init();
	glutDisplayFunc(lineSegment);//ָ����ʾ����
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
