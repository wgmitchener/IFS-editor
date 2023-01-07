IFS-editor
==========

An editor for iterated function system (IFS) fractals.
See *Fractals Everywhere* by Michael Barnsley for the mathematical details of IFS fractals and the multiple-copy process that generates the pictures.

License
-------

This program is free software: you can redistribute it and/or
modify it under the terms of the GNU General Public License as
published by the Free Software Foundation, either version 3 of
the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program, in a file named COPYING.  If not, see
<http://www.gnu.org/licenses/>.


Running the IFS editor
----------------------

The only file you really need is `ifs-editor.rkt`.

To run it, install [Racket](http://www.racket-lang.org).
Launch the DrRacket GUI and load `ifs-editor.rkt`.
Click the run button in the toolbar.

To enter an IFS visually, click on one of the affine transform widgets.
It determines how the unit square, outlined in black, is transformed.
The square handle moves the widget.
The purple circle handle scales and rotates the widget.
The orange circle handle introduces non-uniform scaling and shear.
You can also move the orange circle handle around past the purple one and introduce orientation reversal, which changes the widget's overall color from pink to blue.

The text field in the upper right shows `x0`, `y0`, `a`, `b`, `c`, and `d`, which are the numbers that specify the affine transformation corresponding to currently selected widget.
The transformation it represents is `f(x,y) = A ((x, y) + (x0, y0))` where `A` is the matrix `[a, b; c, d]`.
The `det` line gives the determinant of `A`.
The `fp` line gives the fixed point, that is, the solution to `f(x,y) == (x,y)`.

Click the Edit button under the text field to open a dialog box where you can type in the parameters of the transformation.

The preview area in the lower right shows a couple of iterations of the multiple-copy process that generates the fractal.

Buttons in the toolbar do the following:

- Clear: Delete all affine transform widgets.
- Load: Load a list of transforms from a file.
- Save: Save the current transforms to a file.
- Add: Create a new affine transform widget.
- Delete: Delete the currently selected widget.
- Sketch: Draw a detailed but relatively quick sketch of the invariant set of the IFS.
  (This is the fractal of interest.)
- Render: Draw a detailed rendering of the invariant set of the IFS and save it to an image file.

  - Unit square pixel size: How many pixels correspond to each side of the unit square in the resulting image.
  - Padding pixels: How many more pixels to add on all sides of the unit square to form the resulting image.
  - Maximum recursion depth: How many iterations of the multiple-copy process to perform.

  The resulting image file is square.
  The number of pixels in each dimension is the unit square pixel size plus twice the padding pixels.

For a video demonstration, see <https://youtu.be/kDR9RYSJgjg>.
