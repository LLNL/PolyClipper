// We'll draw a hexagon with one vertex set highlighted

size(150);
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

  // Draw the arrow
  if (i >= vann-1 && i <= vann) {
    draw(a--b, red+linewidth(1.5), MidArrow(arrow_size));
    dot(a, red+linewidth(dot_size));
  } else {
    draw(a--b, black, MidArrow(arrow_size));
    dot(a, black+linewidth(dot_size));
  }

  // Put a labeled point as appropriate.
  real t = pi + (i+1)*dtheta;
  pair labpos = (2*cos(t), 2*sin(t));
  label(format(i), a, labpos);
  if (i == vann) {
    dot(a, blue+linewidth(dot_size));
  } else if (i >= vann-1 && i <= vann+1) {
    dot(a, red+linewidth(dot_size));
  } else {
    dot(a, black+linewidth(dot_size));
  }

  a = b;
  verts.push(b);
}
