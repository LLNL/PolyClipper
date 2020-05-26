// We'll draw a hexagon with one vertex set highlighted

size(400);
import geometry;

real A = 1;       // edge length
int  nfaces = 6;  // hexagon
int  vann = 2;    // Vertex were going to call out and annotate

real dot_size = 10;
real arrow_size = 10;
real line_size = 2;

pair[] verts = {(0,0), (4,0), (4,2), (3,2), (2,1), (1,2), (0,2)};  // Taken from the 2D tests
verts.cyclic = true;                                               // indexing is periodic!

pair[] coord_pos = {NE, NW, SW, SE, S, SW, SE};
pair[] lab_pos = {SW, SE, NE, N, N, N, NW};
for (int i = 0; i < verts.length; ++i) {
  draw(verts[i]--verts[i+1], linewidth(line_size));
  dot(verts[i], black+linewidth(dot_size));
  label(format(i), verts[i], 3*lab_pos[i]);
  label("(" + format(verts[i].x) + ", " + format(verts[i].y) + ")", verts[i], 3*coord_pos[i], magenta);
}

