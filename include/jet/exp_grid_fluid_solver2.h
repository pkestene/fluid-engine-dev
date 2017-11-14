// Copyright (c) 2017 Doyub Kim
//
// I am making my contributions/submissions to this project solely in my
// personal capacity and am not conveying any rights to any intellectual
// property of any third parties.

#ifndef INCLUDE_JET_EXP_GRID_FLUID_SOLVER2_H_
#define INCLUDE_JET_EXP_GRID_FLUID_SOLVER2_H_

#include <jet/array2.h>
#include <jet/physics_animation.h>
#include <jet/size2.h>
#include <jet/vector2.h>

namespace jet {

namespace experimental {

//! Experimental 2-D grid-based fluid solver.
class GridFluidSolver2 : public PhysicsAnimation {
 public:
    //! Default constructor.
    GridFluidSolver2();

    //! Default destructor.
    virtual ~GridFluidSolver2();

    //!
    //! \brief Resizes grid system data.
    //!
    //! This function resizes grid system data. You can also resize the grid by
    //! calling resize function directly from
    //! GridFluidSolver2::gridSystemData(), but this function provides a
    //! shortcut for the same operation.
    //!
    //! \param[in] newSize        The new size.
    //! \param[in] newGridSpacing The new grid spacing.
    //! \param[in] newGridOrigin  The new grid origin.
    //!
    void resizeGrid(const Size2 &newSize, const Vector2D &newGridSpacing,
                    const Vector2D &newGridOrigin);

    ConstArrayAccessor2<float> density() const;

    ArrayAccessor2<float> density();

 protected:
    //! Called when advancing a single time-step.
    void onAdvanceTimeStep(double timeIntervalInSeconds) override;

 private:
    Array2<float> _u;
    Array2<float> _uTemp;
    Array2<float> _v;
    Array2<float> _vTemp;
    Array2<float> _den;
    Array2<float> _denTemp;
};

typedef std::shared_ptr<GridFluidSolver2> GridFluidSolver2Ptr;

}  // namespace experimental

}  // namespace jet

#endif  // INCLUDE_JET_EXP_GRID_FLUID_SOLVER2_H_
