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

class screePt
{
private:
	GLint x;
	GLint y;

public:
	screePt()
	{
		x = y = 0;
	}

	void setCoords(GLint xValue, GLint yValue)
	{
		x = xValue;
		y = yValue;
	}

	GLint getX() const
	{
		return x;
	}

	GLint getY()const
	{
		return y;
	}

	void incrementX()
	{
		x++;
	}

	void decrementY()
	{
		y--;
	}
};

void setPixel(GLint x, GLint y)
{
	glBegin(GL_POINTS);
	glVertex2i(x, y);
	glEnd();
}

void lineSegment(screePt p1, screePt p2)
{
	glBegin(GL_LINES);
	glVertex2i(p1.getX(), p1.getY());
	glVertex2i(p2.getX(), p2.getY());
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

	curvePt[0].setCoords(x0, y0);

	switch (curveNum)
	{
	case limacon:	
		curvePt[0].setCoords(curvePt[0].getX()+ a + b, curvePt[0].getY());
		break;
	case cardioid:	
		curvePt[0].setCoords(curvePt[0].getX() + a + a, curvePt[0].getY());
		break;
	case threeLeaf:	
		curvePt[0].setCoords(curvePt[0].getX() + a, curvePt[0].getY());
		break;
	case fourLeaf:
		curvePt[0].setCoords(curvePt[0].getX() + a, curvePt[0].getY());
		break;
	case spiral: 
		break;
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

		curvePt[1].setCoords(x0 + r * cos(theta), y0 + r * sin(theta));
		lineSegment(curvePt[0], curvePt[1]);

		curvePt[0].setCoords(curvePt[1].getX(), curvePt[1].getY());
		theta += dtheta;
	}

	glEndList();
	return regHex;
}

void initDrawList()
{
	for (int i = 1; i < 6; ++i)
	{
		drawCurve(i);
	}
}

void drawColorTriangles()
{
	glShadeModel(GL_SMOOTH);//GL_SMOOTH-平滑插值，GL_FLAT-单色填充
	glPolygonMode(GL_FRONT, GL_FILL);//GL_LINE-显示边框，GL_POINT-显示顶点，GL_FILL-默认填充显示
	glEnable(GL_LINE_SMOOTH);

	glBegin(GL_TRIANGLES);
	glColor3f(.0f, .0f, 1.f);
	glVertex2i(50, 50);
	glColor3f(.0f, 1.f, .0f);
	glVertex2i(250, 50);
	glColor3f(1.f, .0f, .0f);
	glVertex2i(150, 250);
	glEnd();
}

inline int roundPlus(const float a) { return int(a + 0.5f); }

void lineDDA(int x0, int y0, int xEnd, int yEnd)
{
	int dx = xEnd - x0;
	int dy = yEnd - y0;
	int steps, k;

	float xIncrement, yIncrement, x = x0, y = y0;

	if (abs(dy) > abs(dy))
	{
		steps = dx;
	} 
	else
	{
		steps = dy;
	}

	xIncrement = dx / steps;
	yIncrement = dy / steps;

	setPixel(roundPlus(x), roundPlus(y));
	for (k = 0; k < steps; ++k)
	{
		x = x0 + xIncrement * k;
		y = y0 + yIncrement * k;
		setPixel(roundPlus(x), roundPlus(y));
	}
}

void lineBres(int x0, int y0, int xEnd, int yEnd)
{
	int dx = abs(xEnd - x0);
	int dy = abs(yEnd - y0);
	int dy2 = 2 * dy;
	int a = 2 * (dy - dx);
	int p = dy2 - dx;
	int x, y;

	if (x0 > xEnd)
	{
		x = xEnd;
		y = yEnd;
	} 
	else
	{
		x = x0;
		y = y0;
	}

	setPixel(x, y);
	for (int k = 0; k < dx; ++k)
	{
		x += 1;
		if (p < 0)
		{
			p = p + dy2;
		}
		else
		{
			y += 1;
			p = p + a;
		}
		setPixel(x, y);
	}
}

//中点法画圆
void circlePlotPoints(GLint x, GLint y, screePt point)
{
	setPixel(x + point.getX(), y + point.getY());
	setPixel(x + point.getX(), y - point.getY());
	setPixel(x - point.getX(), y + point.getY());
	setPixel(x - point.getX(), y - point.getY());
	setPixel(x + point.getY(), y + point.getX());
	setPixel(x + point.getY(), y - point.getX());
	setPixel(x - point.getY(), y + point.getX());
	setPixel(x - point.getY(), y - point.getX());
}
void circleMidPoint(GLint x, GLint y, GLint radius)
{
	screePt circPt;
	GLint p = 1 - radius;
	circPt.setCoords(0, radius);
	circlePlotPoints(x, y, circPt);

	while (circPt.getX() < circPt.getY())
	{
		circPt.incrementX();

		if (p < 0)
		{
			p += 2 * circPt.getX() + 1;
		} 
		else
		{
			circPt.decrementY();
			p += 2 * circPt.getX() - 2 * circPt.getY() + 1;
		}

		circlePlotPoints(x, y, circPt);
	}
}

//中点画椭圆
void ellipsePlotPoints(GLint xCenter, GLint yCenter, GLint x, GLint y)
{
	setPixel(xCenter + x, yCenter + y);
	setPixel(xCenter + x, yCenter - y);
	setPixel(xCenter - x, yCenter + y);
	setPixel(xCenter - x, yCenter - y);
}

void ellipseMidPoint(GLint xCenter, GLint yCenter, GLint xRadius, GLint yRadius)
{
	const int Rx2 = xRadius * xRadius;
	const int Ry2 = yRadius * yRadius;
	const int tRx2 = 2 * Rx2;
	const int tRy2 = 2 * Ry2;

	int p = 0;
	int x = 0;
	int y = yRadius;
	int px = 0;
	int py = tRx2 * y;

	ellipsePlotPoints(xCenter, yCenter, x, y);

	p = roundPlus(Ry2 - (Rx2 * yRadius) + (0.25f * Rx2));
	while (px < py)
	{
		x++;
		px += tRy2;
		if (p < 0)
		{
			p += Ry2 + px;
		} 
		else
		{
			y--;
			py -= tRx2;
			p += Ry2 + px - py;
		}
		ellipsePlotPoints(xCenter, yCenter, x, y);
	}

	p = roundPlus(Ry2 * (x + 0.5) * (x + 0.5) + Rx2 * (y - 1) * (y - 1) - Rx2 * Ry2);
	while (y>0)
	{
		y--;
		py -= tRx2;
		if (p > 0)
		{
			p += Rx2 - py;
		}
		else
		{
			x++;
			px += tRy2;
			p += Rx2 - py + px;
		}
		ellipsePlotPoints(xCenter, yCenter, x, y);
	}
}

void midPointGravity(GLint vx0, GLint vy0, GLint x0, GLint y0)
{
	const int gravity = 980;
	const int vxy = vx0 * vy0;
	const int gx = gravity * (x0 - 0.5f);
	const int vxt = vx0 * vx0;
	const int xEnd = (vxy / gravity) + x0;

	int x = x0, y = y0;
	int xg = x * gravity;
	setPixel(x, y);
	setPixel(xEnd + xEnd - x, y);

	int temp = roundPlus((vy0 * vy0 - vx0 * vx0) / (2 * gravity)) + y0;
	int p = vxt + 0.125 * gravity - 0.5 * vxy;
	while (y < temp)
	{
		y++;
		if (p>0)
		{
			x++;
			xg += gravity;
			p += vxt + xg - gx - gravity - vxy;
		} 
		else
		{
			p += vxt;
		}

		setPixel(x, y);
		setPixel(xEnd + xEnd - x, y);
	}
	
	p = (0.5f + y - y0) * vxt - vxy * (x + 1 - x0) + 0.5f * gravity * (x + 1 - x0) * (x + 1 - x0);
	while (x < xEnd)
	{
		x++;
		xg += gravity;
		if (p < 0)
		{
			y++;
			p += xg - gx - vxy + vxt;
		}
		else
		{
			p += xg - gx - vxy;
		}
		setPixel(x, y);
		setPixel(xEnd + xEnd - x, y);
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
	glColor3f(0.0f, 0.4f, 0.2f);//设定显示对象的颜色

	//lineSegment();
	//cube();
	//cubeElements();
	//callList();
	//lineGraph();
	//barChart();
	//drawColorTriangles();
	//lineBres(20, 10, 330, 218);
	//circleMidPoint(200, 200, 100);
	//ellipseMidPoint(200, 200, 100, 50);
	midPointGravity(200, 800, 40, 40);

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
	glViewport(0, 0, newWidth, newHeight);

	glClear(GL_COLOR_BUFFER_BIT);
}

int main(int argc, char** argv) {
	glutInit(&argc, argv);//初始化

	glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);//设定窗口的缓存和颜色模型, GLUT_RGB-颜色编码，GLUT_INDEX-颜色表

	glutInitWindowPosition(100, 100);//设定窗口位置

	glutInitWindowSize(400, 400);//设定窗口大小

	glutCreateWindow("OpenGL Explore");//设置窗口标题

	init();
	initDrawList();
	//glutDisplayFunc(inputDisplayFunc);//指定显示内容
	glutDisplayFunc(drawFunc);//
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
