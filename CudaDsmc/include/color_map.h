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

#ifndef COLOR_MAP_H_
#define COLOR_MAP_H_

#include "config.h"

namespace	dsmc
{
	//create either on gpu or cpu
	//the same to cpu and gpu implementation
	//characterized by strategy

	struct	null_data_provider_t
	{
		void	lock_data(uint width, uint height) {}
		void	unlock_data() {}
		void	init(uint w, uint h) {}
		void	deinit() {}
	};

	struct	gpu_data_provider_t
	{
		void	init(uint w, uint h);
		void	lock_data(uint width, uint height);
		void	unlock_data();
		void	deinit();

		uint	colorBuffer;
	};

	template	<class	DataProvider>
	class	ColorMap
	{
	public:
		void	Deinit();
		void	Render(uint x, uint y, uint w, uint h);
		void	Init(uint width, uint height);
		DataProvider&	GetDataProvider() { return m_dataProvider; }
	private:
		DataProvider	m_dataProvider;
		uint	texId;
		uint	m_width;
		uint	m_height;
	};
}

#include "color_map.inl"

#endif /* COLOR_MAP_H_ */
