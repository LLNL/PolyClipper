// This example shows a cube in PolyClipper topology

import geometry;
import three;

size3(150);

real A = 1;       // edge length
int vann = 4;     // Vertex we're going to call out
int vback = 0;    // Vertex behind the cube

real dot_size = 10;
real arrow_size = 20;

// Pyramid geometry
triple[] coords = {(0,0,0), (A,0,0), (A,A,0), (0,A,0),  // Base
                   (0.5*A, 0.5*A ,A)};                  // Apex
int[][] cube_neighbors = {{1, 4, 3},
                          {2, 4, 0},
                          {3, 4, 1},
                          {0, 4, 2},
                          {0, 1, 2, 3}};

// Draw the edges and such
for (int i = 0; i < coords.length; ++i) {

  // Draw the vertex dots
  bool bluedot = (i == vann);
  bool reddot = false;
  for (int kk = 0; kk < cube_neighbors[vann].length; ++kk) {
    reddot = reddot || (i == cube_neighbors[vann][kk]);
  }
  if (bluedot) {
    dot(coords[i], blue+linewidth(dot_size));
  } else if (reddot) {
    dot(coords[i], red+linewidth(dot_size));
  } else {
    dot(coords[i], linewidth(dot_size));
  }
  label(format(i), coords[i], 2*E);

  // Draw the lines between vertices
  for (int k = 0; k < cube_neighbors[i].length; ++k) {
    int j = cube_neighbors[i][k];

    bool redline = (i == vann || j == vann);
    for (int kk = 0; kk < cube_neighbors[vann].length; ++kk) {
      redline = redline || (j == vann);
    }
    bool dashedline = (i == vback || j == vback);
    
    if (dashedline && redline) {
      draw(coords[i]--coords[j], dashed+red+linewidth(1));
    } else if (dashedline) {
      draw(coords[i]--coords[j], dashed+linewidth(1));
    } else if (redline) {
      draw(coords[i]--coords[j], solid+red+linewidth(1.5));
    } else {
      draw(coords[i]--coords[j], solid+linewidth(1.5));
    }
  }
}

// Draw a directional indicator for how we specify points around a vertex.
triple[] pts = {0.6*coords[vann] + 0.4*coords[cube_neighbors[vann][0]],
                0.6*coords[vann] + 0.4*coords[cube_neighbors[vann][1]],
                0.6*coords[vann] + 0.4*coords[cube_neighbors[vann][2]],
                0.6*coords[vann] + 0.4*coords[cube_neighbors[vann][3]]};
draw(pts[0]..pts[1]..pts[2]..pts[3], Arrow3(arrow_size));

currentprojection=perspective(
camera=(2.96554823701407,5.26167532052144,2.88816390168057),
up=(-0.00156112565675572,-0.00286693866807605,0.00686845868317194),
target=(0.527280331122584,0.531057161148581,0.359381656890007),
zoom=0.872308215245426,
angle=15.1450337404379,
viewportshift=(-0.0161937341305386,0.0685370559208996),
autoadjust=false);
