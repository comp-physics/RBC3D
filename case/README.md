# RBC Example Case

Right now, if you visualize the cells for this case you will see 8 RBC's in a file because that's what the place cells loop in `initoncd.F90` adds.

## Different Cell Types
Currently, the different cell types we have are Leukocytes aka WBC's and sickle cells. We also have functions to create shapes like spheres and ellipsoids in `ModRbc.F90`.

## Adding WBC's
One way to add WBC's:
1. In `initcond.F90`, set rbc%celltype to 2 for the cell you want to make a WBC in and pass that rbc parameter `RBC_Create()` and `RBC_MakeLeukocyte()`.
2. In `initcond.F90`, edit the reference cell array initialization so that rbcRefs at 2 is initialized with `RBC_MakeLeukocyte()`–– this might already be done.

<!--need explanation for why we need ref cell -->

<!-- Why is a function for WBC's called RBC_Create too? Maybe we should rename that -->

## Adding Sickle Cells
