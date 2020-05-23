// This example shows a cube in PolyClipper topology

import geometry;
import three;

size3(150);

real A = 1;       // edge length
int vann = 6;     // Vertex we're going to call out
int vback = 0;   // Vertex behind the cube

real dot_size = 10;
real arrow_size = 10;

// Cube geometry
triple[] coords = {(0,0,0), (A,0,0), (A,A,0), (0,A,0),
                   (0,0,A), (A,0,A), (A,A,A), (0,A,A)};
int[][] cube_neighbors = {{1, 4, 3},
                          {5, 0, 2},
                          {3, 6, 1},
                          {7, 2, 0},
                          {5, 7, 0},
                          {1, 6, 4},
                          {5, 2, 7},
                          {4, 6, 3}};

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
  label(format(i), coords[i], 2*SE);

  // Draw the lines between vertices
  for (int k = 0; k < cube_neighbors[i].length; ++k) {
    int j = cube_neighbors[i][k];

    bool redline = (i == vann || j == vann);
    for (int kk = 0; kk < cube_neighbors[vann].length; ++kk) {
      redline = redline || (j == vann);
    }
    bool dashedline = (i == vback || j == vback);
    
    if (dashedline) {
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
                0.6*coords[vann] + 0.4*coords[cube_neighbors[vann][2]]};
draw(pts[0]..pts[1]..pts[2], Arrow3);
