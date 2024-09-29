import numpy as np
import hist
from hist import Hist


# # Function to rebin a histogram
# def rebin(histogram, new_bins=None, reltol=None, abstol=None):
#     # Check if at least one parameter is provided for rebinning
#     if new_bins is None and reltol is None and abstol is None:
#         raise ValueError("Either 'new_bins', 'reltol', or 'abstol' must be provided for rebinning.")

#     # Process the new bins if provided
#     if new_bins is not None:
#         if not isinstance(new_bins, list):
#             raise ValueError("'new_bins' must be a list or a list of lists.")
#         if not all(isinstance(edges, list) for edges in new_bins):
#             raise ValueError("If 'new_bins' is a list of lists, each element must be a list of bin edges.")
#         # Create a dictionary mapping axis names to the new bin edges
#         new_axes_dict = {axis.name: np.asarray(edges) for axis, edges in zip(histogram.axes, new_bins)}
#     else:
#         new_axes_dict = {}

#     # Process each axis of the histogram
#     for axis in histogram.axes:
#         # If the axis name is not in the new_axes_dict, calculate new edges based on reltol or abstol
#         if axis.name not in new_axes_dict:
#             if reltol is not None and abstol is not None:
#                 raise ValueError("Only one of 'reltol' or 'abstol' should be provided.")
#             elif reltol is not None:
#                 # Calculate new edges based on relative tolerance
#                 rel_edges = np.diff(axis.edges) * reltol
#                 new_edges = axis.edges[:-1][np.where(rel_edges >= 0)]
#                 new_edges = np.concatenate([[axis.edges[0]], new_edges, [axis.edges[-1]]])
#             elif abstol is not None:
#                 # Calculate new edges based on absolute tolerance
#                 new_edges = axis.edges[np.where(np.diff(axis.edges) >= abstol)]
#                 new_edges = np.concatenate([[axis.edges[0]], new_edges, [axis.edges[-1]]])
#             else:
#                 # No rebinning needed for this axis
#                 new_edges = axis.edges
#             # Update the new_axes_dict with the new edges
#             new_axes_dict[axis.name] = new_edges

#     # Calculate the new bin indices for each axis
#     new_axes_indices = {
#         axis.name: np.digitize(histogram.axes[axis.name].centers, new_axes_dict[axis.name]) - 1
#         for axis in histogram.axes
#     }

#     # Define the slices for extracting the rebinned counts
#     slices = tuple(
#         slice(new_axes_indices[axis.name][i], new_axes_indices[axis.name][i+1] + 1)
#         for i, axis in enumerate(histogram.axes)
#     )

#     # Extract the rebinned counts from the original histogram
#     rebinned_counts = histogram.view(flow=True)[slices]

#     # If no counts are present in the rebinned histogram, create an empty array with the correct shape
#     if rebinned_counts.size == 0:
#         rebinned_counts = np.empty(tuple(slice_.stop - slice_.start for slice_ in slices), dtype=histogram.values().dtype)

#     # Create a new rebinned histogram with the updated bin edges
#     rebinned_hist = Hist(
#         *(hist.axis.Variable(new_axes_dict[axis.name]) if axis.name in new_axes_dict else axis for axis in histogram.axes)
#     )

#     return rebinned_hist

# # Function to test the rebin function
# def test_rebin():
#     h = Hist(hist.axis.Regular(5, -5, 5))

#     # Test Case 1
#     rebinned_hist1 = rebin(h, new_bins=[[-7, -5, 0, 1]], reltol=0.1, abstol=None)
#     assert len(rebinned_hist1.axes) == len([[-7, -5, 0, 1]])
#     # Check if the bin edges of the rebinned histogram match the expected values
#     for axis, edges in zip(rebinned_hist1.axes, [[-7, -5, 0, 1]]):
#         assert np.allclose(axis.edges, np.asarray(edges))
#     # Check if the shape of the rebinned histogram values matches the original histogram
#     assert len(rebinned_hist1.values().shape) == len(h.values().shape)
#     assert rebinned_hist1.values().shape == rebinned_hist1.view().shape
#     # Check if the sum of rebinned histogram values matches the sum of original histogram values
#     assert np.allclose(rebinned_hist1.values().sum(), h.values().sum())

#     # Test Case 2
#     rebinned_hist2 = rebin(h, new_bins=[[0, 1, 2.5, 4]], reltol=None, abstol=0.5)
#     assert len(rebinned_hist2.axes) == len([[0, 1, 2.5, 4]])
#     # Check if the bin edges of the rebinned histogram match the expected values
#     for axis, edges in zip(rebinned_hist2.axes, [[0, 1, 2.5, 4]]):
#         assert np.allclose(axis.edges, np.asarray(edges))
#     # Check if the shape of the rebinned histogram values matches the original histogram
#     assert len(rebinned_hist2.values().shape) == len(h.values().shape)
#     assert rebinned_hist2.values().shape == rebinned_hist2.view().shape
#     # Check if the sum of rebinned histogram values matches the sum of original histogram values
#     assert np.allclose(rebinned_hist2.values().sum(), h.values().sum())

#     # Test Case 3
#     rebinned_hist3 = rebin(h, new_bins=[[-1, 0, 1, 2]], reltol=0.2, abstol=0.1)
#     assert len(rebinned_hist3.axes) == len([[-1, 0, 1, 2]])
#     # Check if the bin edges of the rebinned histogram match the expected values
#     for axis, edges in zip(rebinned_hist3.axes, [[-1, 0, 1, 2]]):
#         assert np.allclose(axis.edges, np.asarray(edges))
#     # Check if the shape of the rebinned

def rebin(h, axis_name, edges):
    if type(edges) == int:
        return h[{axis_name : hist.rebin(edges)}]

    ax = h.axes[axis_name]
    ax_idx = [a.name for a in h.axes].index(axis_name)
    if not all([np.isclose(x, ax.edges).any() for x in edges]):
        raise ValueError(f"Cannot rebin histogram due to incompatible edges for axis '{ax.name}'\n"
                            f"Edges of histogram are {ax.edges}, requested rebinning to {edges}")

    # If you rebin to a subset of initial range, keep the overflow and underflow
    overflow = ax.traits.overflow or (edges[-1] < ax.edges[-1] and not np.isclose(edges[-1], ax.edges[-1]))
    underflow = ax.traits.underflow or (edges[0] > ax.edges[0] and not np.isclose(edges[0], ax.edges[0]))
    flow = overflow or underflow
    new_ax = hist.axis.Variable(edges, name=ax.name, overflow=overflow, underflow=underflow)
    axes = list(h.axes)
    axes[ax_idx] = new_ax

    hnew = hist.Hist(*axes, name=h.name, storage=h._storage_type())

    # Offset from bin edge to avoid numeric issues
    offset = 0.5*np.min(ax.edges[1:]-ax.edges[:-1])
    edges_eval = edges+offset
    edge_idx = ax.index(edges_eval)
    # Avoid going outside the range, reduceat will add the last index anyway
    if edge_idx[-1] == ax.size+ax.traits.overflow:
        edge_idx = edge_idx[:-1]

    if underflow:
        # Only if the original axis had an underflow should you offset
        if ax.traits.underflow:
            edge_idx += 1
        edge_idx = np.insert(edge_idx, 0, 0)

    # Take is used because reduceat sums i:len(array) for the last entry, in the case
    # where the final bin isn't the same between the initial and rebinned histogram, you
    # want to drop this value. Add tolerance of 1/2 min bin width to avoid numeric issues
    hnew.values(flow=flow)[...] = np.add.reduceat(h.values(flow=flow), edge_idx,
            axis=ax_idx).take(indices=range(new_ax.size+underflow+overflow), axis=ax_idx)
    if hnew._storage_type() == hist.storage.Weight():
        hnew.variances(flow=flow)[...] = np.add.reduceat(h.variances(flow=flow), edge_idx,
                axis=ax_idx).take(indices=range(new_ax.size+underflow+overflow), axis=ax_idx)
    return hnew