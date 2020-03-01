// openGL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <Windows.h>
//#include <GL/glew.h>
#include <GL/glut.h>

#include <iostream>

using namespace std;

//#pragma comment(lib, "glew32.lib")

const double TWO_PI = 6.2831853;
GLsizei winWidth = 600, winHeight = 500;
GLint xRaster = 25, yRaster = 100;
GLubyte label[36] = { 'J', 'a', 'n',		'F', 'e', 'b',		'M', 'a', 'r',
					  'A', 'p', 'r',		'M', 'a', 'y',		'J', 'u', 'n',
					  'J', 'u', 'l',		'A', 'u', 'g',		'S', 'e', 'p',
					  'O', 'c', 't',		'N', 'o', 'v',		'D', 'e', 'c',
};

GLint dataValue[12] = { 420, 342, 324, 310, 262, 185, 190, 196, 217, 240, 312, 438 };

struct screePt
{
	GLint x;
	GLint y;
};

void lineSegment(screePt p1, screePt p2)
{
	glBegin(GL_LINES);
	glVertex2i(p1.x, p1.y);
	glVertex2i(p2.x, p2.y);
	glEnd();
}

typedef GLint vertex3[3];
vertex3 pt[8] = { {0, 0, 0}, {0, 100, 0}, {100, 0, 0}, {100, 100, 0}, {0, 0, 100}, {0, 100, 100}, {100, 0, 100}, {100, 100, 100} };
void quad(GLint n1, GLint n2, GLint n3, GLint n4)
{
	glBegin(GL_QUADS);
	glVertex3iv(pt[n1]);
	glVertex3iv(pt[n2]);
	glVertex3iv(pt[n3]);
	glVertex3iv(pt[n4]);
	glEnd();
}

void cube()
{
	quad(6, 2, 3, 7);
	quad(5, 1, 0, 4);
	quad(7, 3, 1, 5);
	quad(4, 0, 2, 6);
	quad(2, 0, 1, 3);
	quad(4, 6, 7, 5);
}

void cubeElements()
{
	glEnableClientState(GL_VERTEX_ARRAY);//激活顶点数组
	glVertexPointer(3, GL_INT, 0, pt);//设置顶点坐标的位置和格式

	GLubyte vertexIndex[] = { 6, 2, 3, 7, 5, 1, 0, 4, 7, 3, 1, 5, 4, 0, 2, 6, 2, 0, 1, 3, 7, 5, 4, 6 };//设置顶点索引
	glDrawElements(GL_QUADS, 24, GL_UNSIGNED_BYTE, vertexIndex);
}

void callList()
{
	GLuint regHex;
	GLdouble theta;
	GLint x, y, k;

	regHex = glGenLists(1);
	glNewList(regHex, GL_COMPILE);
	glBegin(GL_POLYGON);
	for (k = 0; k < 6; ++k)
	{
		theta = TWO_PI * k / 6.0;
		x = 200 + 150 * cos(theta);
		y = 200 + 150 * sin(theta);
		glVertex2i(x, y);
	}
	glEnd();
	glEndList();

	glCallList(regHex);
}

void lineGraph()
{
	GLint month, k;
	GLint x = 30;

	glColor3f(0.0f, 0.0f, 1.0f);
	glBegin(GL_LINE_STRIP);
	for (k = 0; k <12;k++)
	{
		glVertex2i(x + k * 50, dataValue[k]);
	}
	glEnd();

	glColor3f(1.0f, 0.0f, 0.0f);
	xRaster = 25;
	for (k = 0; k < 12; k++)
	{
		glRasterPos2i(xRaster + k * 50, dataValue[k] - 4);
		glutBitmapCharacter(GLUT_BITMAP_9_BY_15, '*');
	}

	glColor3f(0.0f, 0.0f, 0.0f);
	xRaster = 20;
	for (month = 0; month < 12; month++)
	{
		glRasterPos2i(xRaster, yRaster);
		for (k = 3*month; k < 3*month+ 3; k++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[k]);
		}
		xRaster += 50;
	}
}

void barChart()
{
	GLint month, k;

	glColor3f(1.0, 0.0, 0.0);
	for (k = 0; k < 12; ++k)
	{
		glRecti(20 + 50 * k, 150, 40 + 50 * k, dataValue[k]);
	}

	glColor3f(0.0f, 0.0f, 0.0f);
	xRaster = 20;
	for (month = 0; month < 12; month++)
	{
		glRasterPos2i(xRaster, yRaster);
		for (k = 3 * month; k < 3 * month + 3; k++)
		{
			glutBitmapCharacter(GLUT_BITMAP_HELVETICA_12, label[k]);
		}
		xRaster += 50;
	}
}

typedef enum 
{
	limacon = 1,
	cardioid,
	threeLeaf,
	fourLeaf,
	spiral,
} curveNames;

GLuint drawCurve(GLint curveNum)
{
	GLuint regHex = glGenLists(1);
	glNewList(regHex, GL_COMPILE);

	const GLint a = 175, b = 60;

	GLfloat r, theta, dtheta = 1.0 / (float)(a);
	GLint x0 = 200, y0 = 250;

	screePt curvePt[2];

	glColor3f(0.0f, 0.0f, 0.0f);

	curvePt[0].x = x0;
	curvePt[0].y = y0;

	switch (curveNum)
	{
	case limacon:	curvePt[0].x += a + b; break;
	case cardioid:	curvePt[0].x += a + a; break;
	case threeLeaf:	curvePt[0].x += a; break;
	case fourLeaf:	curvePt[0].x += a; break;
	case spiral: break;
	default:
		break;
	}

	theta = dtheta;
	while (theta < TWO_PI)
	{
		switch (curveNum)
		{
		case limacon:
			r = a * cos(theta) + b; break;
		case cardioid:
			r = a * (1 + cos(theta)); break;
		case threeLeaf:
			r = a * cos(3 * theta); break;
		case fourLeaf:
			r = a * cos(2 * theta); break;
		case spiral:
			r = (a / 4.0) * theta; break;
		default:
			break;
		}

		curvePt[1].x = x0 + r * cos(theta);
		curvePt[1].y = y0 + r * sin(theta);
		lineSegment(curvePt[0], curvePt[1]);

		curvePt[0].x = curvePt[1].x;
		curvePt[0].y = curvePt[1].y;
		theta += dtheta;
	}

	glEndList();
	return regHex;
}

void initDrawList()
{
	for (int i = 1; i < 6; ++i)
	{
		GLuint regHex = drawCurve(i);

	}
}
void inputDisplayFunc(void)
{
	GLuint curveNum;

	glClear(GL_COLOR_BUFFER_BIT);
	cout << "\n Enter the integer value corresponding to \n";
	cout << "one of the following curve names.\n";
	cout << "Press any other key to exit.\n";
	cout << "\n 1.limacon, 2-cardioid, 3-threeLeaf, 4-fourLeaf, 5-spiral:";
	cin >> curveNum;

	if (glIsList(curveNum) == GL_TRUE)
		glCallList(curveNum);
	else
		exit(0);

	glFlush();
}

void drawFunc(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glColor3f(0.6f, 0.4f, 0.2f);//设定显示对象的颜色

	//lineSegment();
	//cube();
	//cubeElements();
	//callList();
	//lineGraph();
	barChart();

	glFlush();
}

void init(void)
{
	glClearColor(1.0f, 1.0f, 1.0f, 0.0f);//设定背景颜色

	glMatrixMode(GL_PROJECTION);//设定投影模式
	glLoadIdentity();
	gluOrtho2D(0.0, 400.0, 0, 400.0);
}

void winReshapeFunc(int newWidth, int newHeight)
{
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluOrtho2D(0, (GLdouble)newWidth, 0, (GLdouble)newHeight);
	//glViewport(0, 0, newWidth, newHeight);

	glClear(GL_COLOR_BUFFER_BIT);
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);//初始化

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);//设定窗口的缓存和颜色模型

	glutInitWindowPosition(100, 100);//设定窗口位置

	glutInitWindowSize(400, 400);//设定窗口大小

	glutCreateWindow("OpenGL Explore");//设置窗口标题

	init();
	initDrawList();
	glutDisplayFunc(inputDisplayFunc);//指定显示内容
	glutReshapeFunc(winReshapeFunc);//显示窗口重定形
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
