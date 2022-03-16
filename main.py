import pyrvgomea
x0, y0, x1, y1 = 6, 2, 3, 4
rect_obj = pyrvgomea.pyrvgomea(x0, y0, x1, y1)
print(dir(rect_obj))
print(rect_obj.get_size())
print(rect_obj.get_area())
