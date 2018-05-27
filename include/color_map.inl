/*
The MIT License (MIT)

Copyright (c) 2013 Denis Gladkov

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#include "draw_helper.h"

#include <GL/glew.h>

namespace	dsmc
{

template	<class	DataProvider>
void	ColorMap<DataProvider>::Render(uint x, uint y, uint w, uint h)
{
	glBindTexture(GL_TEXTURE_2D,texId);

	glColor3f(1,1,1);

	m_dataProvider.lock_data(m_width, m_height);

	glBegin(GL_QUADS);

		glTexCoord2f(0,0);
		glVertex2i(x,y);

		glTexCoord2f(0,1);
		glVertex2i(x,y+h);

		glTexCoord2f(1,1);
		glVertex2i(x+w, y+h);

		glTexCoord2f(1,0);
		glVertex2i(x+w, y);

	glEnd();

	m_dataProvider.unlock_data();

	glBindTexture(GL_TEXTURE_2D, 0);

}

template	<class	DataProvider>
void	ColorMap<DataProvider>::Deinit()
{
	m_dataProvider.deinit();
}

template	<class	DataProvider>
void	ColorMap<DataProvider>::Init(uint width, uint height)
{
	m_width = width;
	m_height = height;

	texId = CreateTexture(m_width, m_height);
	m_dataProvider.init(m_width, m_height);
}

/*
		DataProvider	m_dataProvider;
		uint	texId;
*/
}
