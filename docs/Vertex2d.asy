// We'll draw a hexagon with one vertex set highlighted

size(200);
import geometry;

real A = 1;       // edge length
int  nfaces = 6;  // hexagon
int  vann = 2;    // Vertex were going to call out and annotate

real dot_size = 10;
real arrow_size = 10;
real dtheta = 2.0*pi/nfaces;
pair a = (0, 0);
pair[] verts = {a};

for (int i = 0; i < nfaces; ++i) {
  real theta = i*dtheta;
  pair b = a + (A*cos(theta), A*sin(theta));
  string lab = format(i);
  if (i >= vann-1 && i <= vann) {
    draw(Label(lab, BeginPoint, black), a--b, red+linewidth(1.5), MidArrow(arrow_size));
    dot(a, red+linewidth(dot_size));
  } else {
    draw(Label(lab, BeginPoint, black), a--b, black, MidArrow(arrow_size));
    dot(a, black+linewidth(dot_size));
  }
  a = b;
  verts.push(b);
}

// dot(verts[vann-1], red+linewidth(dot_size));
// dot(verts[vann],   blue+linewidth(dot_size));
// dot(verts[vann+1], red+linewidth(dot_size));
