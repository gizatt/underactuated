# -*- coding: utf8 -*-

import argparse
import math
import os.path

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy as sp

from drake import lcmt_viewer_load_robot
from pydrake.common.eigen_geometry import Quaternion
from pydrake.geometry import DispatchLoadMessage, SceneGraph
from pydrake.lcm import DrakeMockLcm
from pydrake.math import RigidTransform, RotationMatrix
from pydrake.systems.rendering import PoseBundle

from pydrake.all import (
    AbstractValue,
    AddModelInstanceFromUrdfStringSearchingInRosPackages,
    ConstantVectorSource,
    Context,
    DiagramBuilder,
    FloatingBaseType,
    Parser,
    PortDataType,
    MultibodyPlant,
    SceneGraph,
    SignalLogger,
    Simulator,
)

from utils import FindResource
from pyplot_visualizer import PyPlotVisualizer


class PlanarMultibodyVisualizer(PyPlotVisualizer):

    '''
    Given a SceneGraph and a view plane, provides a view of the robot by
    projecting all geometry onto the view plane.

    This is intended to be used for robots that operate in the plane, and won't
    render out-of-plane transformations correctly. (It precomputes how geometry
    looks on the projected plane and assumes it won't be rotated out of plane.
    It'll *work* for out-of-plane transformations, but render things
    increasingly inaccurately.)

    Params:
    - Tview, xlim, and ylim set up view into scene.
    - facecolor is passed through to figure() and sets background color. Both
    color name strings and RGB triplets are allowed. Defaults to white.
    - use_random_colors, if set to True, will render each body with a different
    color. (Multiple visual elements on the same body will be the same color.)
    - if ax is supplied, the visualizer will draw onto those axes instead of
    creating a new set of axes. The visualizer will still change the view range
    and figure size of those axes.

     Specifics on view setup:

    TView specifies the view projection matrix,
    and should be a 3x4 matrix:
    [ <x axis select> x_axis_shift
      <y axis select> y_axis_shift
       0, 0, 0, 1]  % homogenizer

    e.g.

    [ 1 0 0 0.5
      0 1 0 0
      0 0 0 1]

    would give a top-down view (i.e squashing the z axis), and would shift
    things in the x axis positively by 0.5.

    xlim and ylim don't technically provide extra functionality, but I think
    it's easier to keep handle scaling with xlim and ylim and view plane
    selection and *maybe* offsetting with the projection matrix.
    '''

    def __init__(self,
                 scene_graph,
                 draw_period=0.033333,
                 Tview=np.array([[1., 0., 0., 0.],
                                 [0., 0., 1., 0.],
                                 [0., 0., 0., 1.]]),
                 xlim=[-1., 1],
                 ylim=[-1, 1],
                 facecolor=[1, 1, 1],
                 use_random_colors=False,
                 ax=None):

        default_size = matplotlib.rcParams['figure.figsize']
        scalefactor = (ylim[1]-ylim[0])/(xlim[1]-xlim[0])
        figsize = (default_size[0], default_size[0]*scalefactor)

        PyPlotVisualizer.__init__(self, facecolor=facecolor, figsize=figsize,
                                  ax=ax, draw_timestep=draw_period)
        self.set_name('planar_multibody_visualizer')

        self._scene_graph = scene_graph
        self.Tview = Tview
        self.Tview_pinv = np.linalg.pinv(self.Tview)

        # Pose bundle (from SceneGraph) input port.
        self._DeclareAbstractInputPort("lcm_visualization",
                                       AbstractValue.Make(PoseBundle(0)))

        self.ax.axis('equal')
        self.ax.axis('off')

        # Achieve the desired view limits
        self.ax.set_xlim(xlim)
        self.ax.set_ylim(ylim)
        default_size = self.fig.get_size_inches()
        scalefactor = (ylim[1]-ylim[0])/(xlim[1]-xlim[0])
        self.fig.set_size_inches(default_size[0],
                                 default_size[0]*scalefactor)

        # Populate body patches
        self.buildViewPatches(use_random_colors)

        # Populate the body fill list -- which requires doing most of
        # a draw pass, but with an ax.fill() command rather than
        # an in-place replacement of vertex positions to initialize
        # the draw patches.
        # The body fill list stores the ax patch objects in the
        # order they were spawned (i.e. by body, and then by
        # order of viewPatches). Drawing the tree should update them
        # by iterating over bodies and patches in the same order.
        self.body_fill_dict = {}
        n_bodies = len(self.viewPatches.keys())
        tf = np.eye(4)
        for full_name in self.viewPatches.keys():
            viewPatches, viewColors = self.getViewPatches(full_name, tf)
            self.body_fill_dict[full_name] = []
            for patch, color in zip(viewPatches, viewColors):
                patch_proj = np.dot(self.Tview, patch)
                self.body_fill_dict[full_name] += self.ax.fill(
                    patch_proj[0, :], patch_proj[1, :], zorder=0, edgecolor='k',
                    facecolor=color, closed=True)

    def buildViewPatches(self, use_random_colors):
        ''' Generates view patches. self.viewPatches stores a list of
        viewPatches for each body (starting at body id 1). A viewPatch is a
        list of 2D coordinates in counterclockwise order forming the boundary
        of a filled polygon representing a piece of visual geometry. '''

        self.viewPatches = {}
        self.viewPatchColors = {}

        mock_lcm = DrakeMockLcm()
        DispatchLoadMessage(self._scene_graph, mock_lcm)
        load_robot_msg = lcmt_viewer_load_robot.decode(
            mock_lcm.get_last_published_message("DRAKE_VIEWER_LOAD_ROBOT"))

        # Spawn a random color generator, in case we need to pick
        # random colors for some bodies. Each body will be given
        # a unique color when using this random generator, with
        # each visual element of the body colored the same.
        color = iter(plt.cm.rainbow(np.linspace(0, 1, load_robot_msg.num_links)))

        for i in range(load_robot_msg.num_links):
            link = load_robot_msg.link[i]

            this_body_patches = []
            this_body_colors = []
            this_color = next(color)

            for j in range(link.num_geom):
                geom = link.geom[j]
                # MultibodyPlant currently sets alpha=0 to make collision
                # geometry "invisible".  Ignore those geometries here.
                if geom.color[3] == 0:
                    continue

                element_local_tf = RigidTransform(
                    RotationMatrix(Quaternion(geom.quaternion)),
                    geom.position)

                if geom.type == geom.BOX:
                    assert geom.num_float_data == 3

                    # Draw a bounding box.
                    patch = np.vstack((
                        geom.float_data[0]/2.*np.array([-1, -1, 1, 1, -1, -1, 1, 1]),
                        geom.float_data[1]/2.*np.array([-1, 1, -1, 1, -1, 1, -1, 1]),
                        geom.float_data[2]/2.*np.array([-1, -1, -1, -1, 1, 1, 1, 1])))

                elif geom.type == geom.SPHERE:
                    assert geom.num_float_data == 1
                    radius = geom.float_data[0]
                    lati, longi = np.meshgrid(np.arange(0., 2.*math.pi, 0.5),
                                              np.arange(0., 2.*math.pi, 0.5))
                    lati = lati.ravel()
                    longi = longi.ravel()
                    patch = np.vstack([
                        np.sin(longi)*np.cos(lati),
                        np.sin(longi)*np.sin(lati),
                        np.cos(lati)])
                    patch *= radius

                elif geom.type == geom.CYLINDER:
                    assert geom.num_float_data == 2
                    radius = geom.float_data[0]
                    length = geom.float_data[1]

                    # In the lcm geometry, cylinders are along +z

                    # https://github.com/RobotLocomotion/drake/blob/last_sha_with_original_matlab/drake/matlab/systems/plants/RigidBodyCylinder.m

                    # Two circles: one at bottom, one at top.
                    sample_pts = np.arange(0., 2.*math.pi, 0.25)
                    patch = np.hstack(
                        [np.array([
                            [radius*math.cos(pt),
                            radius*math.sin(pt),
                             -length/2.],
                            [radius*math.cos(pt),
                             radius*math.sin(pt),
                             length/2.]]).T
                         for pt in sample_pts])

                elif geom.type == geom.MESH:
                    import trimesh
                    mesh = trimesh.load(geom.string_data)
                    patch = mesh.vertices.T
                    # Apply scaling
                    for i in range(3):
                        patch[i, :] *= geom.float_data[i]
                else:
                    print("UNSUPPORTED GEOMETRY TYPE {} IGNORED".format(
                        geom.type))
                    continue

                patch = np.vstack((patch, np.ones((1, patch.shape[1]))))
                patch = np.dot(element_local_tf.GetAsMatrix4(), patch)
                # Project into 2D
                # patch = np.dot(self.Tview, patch)

                # Close path if not closed
                if (patch[:, -1] != patch[:, 0]).any():
                    patch = np.hstack((patch, patch[:, 0][np.newaxis].T))

                this_body_patches.append(patch)
                if use_random_colors:
                    this_body_colors.append(this_color)
                else:
                    this_body_colors.append(geom.color)

            self.viewPatches[link.name] = this_body_patches
            self.viewPatchColors[link.name] = this_body_colors

    def getViewPatches(self, full_name, tf):
        ''' Pulls out the view patch verts for the given body index after
            applying the appropriate TF '''
        #projected_tf = np.dot(np.dot(self.Tview, tf), self.Tview_pinv)
        transformed_patches = [np.dot(tf, patch)
                               for patch in self.viewPatches[full_name]]
        colors = self.viewPatchColors[full_name]
        return (transformed_patches, colors)

    def draw(self, context):
        ''' Evaluates the robot state and draws it.
            Can be passed either a raw state vector, or
            an input context.'''

        pose_bundle = self.EvalAbstractInput(context, 0).get_value()

        for frame_i in range(pose_bundle.get_num_poses()):
            # SceneGraph currently sets the name in PoseBundle as
            #    "get_source_name::frame_name".
            full_name = pose_bundle.get_name(frame_i)
            model_id = pose_bundle.get_model_instance_id(frame_i)

            viewPatches, _ = self.getViewPatches(
                full_name,
                pose_bundle.get_pose(frame_i).matrix())
            for i, patch in enumerate(viewPatches):
                patch_proj = np.dot(self.Tview, patch)[0:2, :]
                if patch_proj.shape[1] > 3:
                    hull = sp.spatial.ConvexHull(
                        np.transpose(patch_proj[0:2, :]))
                    patch_proj = np.transpose(
                        np.vstack([patch_proj[:, v] for v in hull.vertices]))
                # patch_proj[:, -1] = patch_proj[:, 0]
                n_verts = self.body_fill_dict[full_name][i].get_path().vertices.shape[0]
                patch_proj = np.pad(patch_proj, ((0, 0), (0, n_verts - patch_proj.shape[1])), mode="edge")
                self.body_fill_dict[full_name][i].get_path().vertices[:, :] = np.transpose(patch_proj)  # noqa

        if isinstance(context, Context):
            self.ax.set_title('t = {:.1f}'.format(context.get_time()))
        else:
            self.ax.set_title('')


def setupPendulumExample(plant):
    return rbt, pbrv


def setupDoublePendulumExample(plant):
    rbt = RigidBodyTree(FindResource("double_pendulum/double_pendulum.urdf"),
                        floating_base_type=FloatingBaseType.kFixed)  # noqa
    Tview = np.array([[1., 0., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., 1.]],
                     dtype=np.float64)
    fig, ax = plt.subplots(1, 1)
    pbrv = PlanarRigidBodyVisualizer(rbt, Tview, [-2.5, 2.5], [-2.5, 2.5],
                                     use_random_colors=True, ax=ax)
    return rbt, pbrv


def setupValkyrieExample(plant):
    # Valkyrie Example
    rbt = RigidBodyTree()
    world_frame = RigidBodyFrame("world_frame", rbt.world(),
                                 [0, 0, 0], [0, 0, 0])
    from pydrake.multibody.parsers import PackageMap
    import pydrake
    pmap = PackageMap()
    # Note: Val model is currently not installed in drake binary distribution.
    pmap.PopulateFromFolder(os.path.join(pydrake.getDrakePath(), "examples"))
    # TODO(russt): remove plane.urdf and call AddFlatTerrainTOWorld instead
    AddModelInstanceFromUrdfStringSearchingInRosPackages(
        open(FindResource(os.path.join("underactuated", "plane.urdf")), 'r').read(),  # noqa
        pmap,
        pydrake.getDrakePath() + "/examples/",
        FloatingBaseType.kFixed,
        world_frame,
        rbt)
    val_start_frame = RigidBodyFrame("val_start_frame", rbt.world(),
                                     [0, 0, 1.5], [0, 0, 0])
    AddModelInstanceFromUrdfStringSearchingInRosPackages(
        open(pydrake.getDrakePath() + "/examples/valkyrie/urdf/urdf/valkyrie_A_sim_drake_one_neck_dof_wide_ankle_rom.urdf", 'r').read(),  # noqa
        pmap,
        pydrake.getDrakePath() + "/examples/",
        FloatingBaseType.kRollPitchYaw,
        val_start_frame,
        rbt)
    Tview = np.array([[1., 0., 0., 0.],
                      [0., 0., 1., 0.],
                      [0., 0., 0., 1.]],
                     dtype=np.float64)
    pbrv = PlanarRigidBodyVisualizer(rbt, Tview, [-2.0, 2.0], [-0.25, 3.0],
                                     use_random_colors=True)
    return rbt, pbrv


if __name__ == "__main__":
    # Usage demo: load a URDF, rig it up with a constant torque input, and
    # draw it.

    np.set_printoptions(precision=5, suppress=True)
    parser = argparse.ArgumentParser()
    parser.add_argument("-T", "--duration",
                        type=float,
                        help="Duration to run sim.",
                        default=1.0)
    parser.add_argument("-m", "--models",
                        type=str,
                        nargs="*",
                        help="Models to run, at least one of [pend, dpend, "
                             "val]",
                        default=["pend"])
    parser.add_argument("-a", "--animate",
                        action="store_true",
                        help="Enable real-time looping animation after each "
                             "simulation.")
    args = parser.parse_args()

    for model in args.models:
        builder = DiagramBuilder()
        timestep = 0.001
        plant = builder.AddSystem(MultibodyPlant(timestep))
        scene_graph = builder.AddSystem(SceneGraph())
        plant.RegisterAsSourceForSceneGraph(scene_graph)
        builder.Connect(plant.get_geometry_poses_output_port(),
                        scene_graph.get_source_pose_port(
                            plant.get_source_id()))
        parser = Parser(plant)

        if model == "pend":
            parser.AddModelFromFile(FindResource("pendulum/pendulum.urdf"))
            Tview = np.array([[1., 0., 0., 0.],
                              [0., 0., 1., 0.],
                              [0., 0., 0., 1.]],
                             dtype=np.float64)
            visualizer = builder.AddSystem(PlanarMultibodyVisualizer(
                scene_graph, Tview, [-1.2, 1.2], [-1.2, 1.2]))
        elif model == "dpend":
            visualizer = builder.AddSystem(setupDoublePendulumExample(plant))
        elif model == "val":
            visualizer = builder.AddSystem(setupValkyrieExample(plant))
        else:
            print "Unrecognized model %s." % model
            parser.print_usage()
            exit(1)

        plant.Finalize()
        builder.Connect(scene_graph.get_pose_bundle_output_port(),
                        visualizer.get_input_port(0))

        diagram = builder.Build()
        simulator = Simulator(diagram)
        simulator.Initialize()
        simulator.set_target_realtime_rate(1.0)

        # Fix the input port to zero.
        plant_context = diagram.GetMutableSubSystemContext(
            simulator.get_mutable_context(), plant)
        plant_context.FixInputPort(0, np.zeros(plant.num_actuators()))

        simulator.StepTo(args.duration)

        # Generate an animation of whatever happened
        ani = visualizer.animate(signalLogger, repeat=True)

