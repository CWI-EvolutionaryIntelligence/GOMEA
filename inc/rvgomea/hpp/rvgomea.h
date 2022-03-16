#ifndef RECTANGLE_H
#define RECTANGLE_H

namespace shapes {
    class rvgomea {
        public:
            int x0, y0, x1, y1;
            rvgomea();
            rvgomea(int x0, int y0, int x1, int y1);
            ~rvgomea();
            int getArea();
            void getSize(int* width, int* height);
            void move(int dx, int dy);
    };
}

#endif
