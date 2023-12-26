# RBC Example Case

Right now, if you visualize the cells for this case, you will see 8 RBC's in a file because that's what the place cells loop in `initoncd.F90` adds.

## Different Cell Types
As of now, the different cell types we have are Leukocytes (aka WBC's) and sickle cells. We also have functions to create spheres and ellipsoids in `ModRbc.F90`.  

When you place a cell, you have to make sure the cell has a celltype number corresponding to the correct reference cell. For example, a RBC with celltype 1 has to have the reference cell at `rbcRefs(1)` created via `RBC_MakeBiConcave`. This is already done.

## Adding WBC's
One way to make sure WBC's are added correctly:
1. In `initcond.F90`, set rbc%celltype to 2 for the cell you want to make a WBC in and pass that rbc object into `RBC_Create()` and `RBC_MakeLeukocyte()`.
2. In `tube.F90`, edit the reference cell array initialization so that the `rbcRefs(2)` is initialized with `RBC_MakeLeukocyte()`.

## Adding Sickle Cells
