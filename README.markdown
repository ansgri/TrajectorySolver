TrajectorySolver
================

TrajectorySolver is a small Java application with [Cliché](http://code.google.com/p/cliche/)-based UI for computation of the motion of the charged particle in given electrostatic field.

For diagnostic purposes it can draw charts of particle's energy and coordinates (using [JChart2D](http://jchart2d.sourceforge.net/) charting library).

The project aims to provide a simple tool for accurate *energy-conserving* solution of Newton's equations in electrostatic field, though this goal is not yet reached.

Unit system
-----------

<table>
  <tr><th>quantity</th><th>units</th></tr>
  <tr><td>length</td><td>k * meter</td></tr>
  <tr><td>time</td><td>k * second</td></tr>
  <tr><td>energy</td><td>electron-volt</td></tr>
  <tr><td>potential</td><td>volt</td></tr>
  <tr><td>charge</td><td>units of |e|</td></tr>
</table>

