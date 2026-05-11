#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""Plot graph's connected (connex) components.
Takes as input a graph and plots the N biggest connected components.
"""


import argparse
import os
import pickle
import sys
import time

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np

import utils


def _rank_clusters(graph):
    """Sort connected components by total reads, descending."""
    ranked = []
    for component in nx.connected_components(graph):
        sub = graph.subgraph(component)
        total = sum(sub.nodes[n]['n_reads'] for n in sub)
        ranked.append((total, sub))
    ranked.sort(key=lambda t: t[0], reverse=True)
    return ranked


def _draw_cluster(sub, ax, fig, *, node_size, node_label_font,
                  edge_label_font, title_font, node_label_max_nodes,
                  edge_label_max_edges, cmap):
    """Draw a single cluster (connected component) on the given axis.
    Heavier nodes are drawn last so they end up painted on top of lighter
    ones — useful for spotting mass centers in dense graphs.

    Arguments:
    sub             (nx.Graph) - Connected (sub-)graph.
    ax              (plt.Axes) - Matplotlib subplot.
    fig           (plt.Figure) - Matplotlib figure.
    node_size          (float) - Fixed size of nodes in the graph.
    node_label_font    (float) - Size of nodes labels font.
    edge_label_font    (float) - Size of edges labels font.
    title_font         (float) - Size of title.
    node_label_max_nodes (int) - Maximum number of nodes for nodes labels to be
                                 shown.
    edge_label_max_edges (int) - Maximum number of edges for edges labels to be
                                 shown.
    cmap        (plt.Colormap) - Color map used for nodes coloring according to
                                 number of node's reads.

    Returns:
    rep         (str) - Label of representative node i.e. node with greatest
                        weight.
    total_reads (int) - Total amount of reads in the connected (sub-)graph.
    """
    pos = nx.spring_layout(sub, weight=None, seed=0)

    node_reads = {n: sub.nodes[n]['n_reads'] for n in sub}
    # Get extreme values for color map
    vmin = min(node_reads.values())
    vmax = max(node_reads.values())
    if vmax == vmin:
        vmax = vmin + 1  # avoid degenerate single-value colormap

    # Draw edges first (under nodes)
    nx.draw_networkx_edges(sub, pos, ax=ax, alpha=0.45, width=0.8)

    # Split nodes by weight class and draw in two passes. Singletons go in a
    # translucent, smaller back layer; heavier nodes go in an opaque,
    # full-size front layer sorted ascending so the heaviest are painted on
    # top. Splitting by alpha class instead of passing a list of per-node
    # alphas avoids networkx pre-multiplying colors through the cmap (which
    # yields RGBA data and triggers a "no data for colormapping" warning
    # from matplotlib's scatter).
    singletons = [n for n in sub if node_reads[n] == 1]
    heavies = sorted(
        [n for n in sub if node_reads[n] > 1],
        key=lambda n: node_reads[n],
    )

    # Back layer: translucent singletons (skipped if there are no heavies to
    # contrast with -- in that case singletons are drawn opaque below)
    if singletons and heavies:
        nx.draw_networkx_nodes(
            sub, pos, ax=ax,
            nodelist=singletons,
            node_size=node_size * 0.4,
            node_color=[node_reads[n] for n in singletons],
            alpha=0.35,
            cmap=cmap, vmin=vmin, vmax=vmax,
            edgecolors='#333333', linewidths=0.5,
        )

    # Front layer: heavies on top, sorted ascending. Falls back to drawing
    # the singletons here if the cluster has no heavy node, so the colorbar
    # still has a valid mappable.
    front = heavies if heavies else singletons
    nodes = nx.draw_networkx_nodes(
        sub, pos, ax=ax,
        nodelist=front,
        node_size=node_size,
        node_color=[node_reads[n] for n in front],
        alpha=1.0,
        cmap=cmap, vmin=vmin, vmax=vmax,
        edgecolors='#333333', linewidths=0.5,
    )

    # Draw nodes labels
    if sub.number_of_nodes() <= node_label_max_nodes:
        nx.draw_networkx_labels(
            sub, pos, ax=ax, font_size=node_label_font,
            labels={n: node_reads[n] for n in sub},
        )
    
    # Old way of scaling node sizes with number of reads
    # sizes = [80 + 25 * sub.nodes[n]['n_reads'] for n in sub]
    # nx.draw_networkx_nodes(sub, pos, ax=ax, node_size=sizes,
    #                        node_color='#4C9CD6', alpha=0.85)
    # nx.draw_networkx_edges(sub, pos, ax=ax, alpha=0.5)
    # nx.draw_networkx_labels(
    #     sub, pos, ax=ax, font_size=7,
    #     labels={n: sub.nodes[n]['n_reads'] for n in sub},
    # )
    
    # Draw edges labels
    if sub.number_of_edges() <= edge_label_max_edges:
        edge_labels = {
            (u, v): f"{d['alignment_distance']}|"
                    f"{d['overlap_distance']:.2f}"
            for u, v, d in sub.edges(data=True)
        }
        nx.draw_networkx_edge_labels(
            sub, pos, edge_labels=edge_labels, ax=ax,
            font_size=edge_label_font,
        )

    # Get node with maximum amount of reads
    rep = max(sub.nodes, key=lambda n: node_reads[n])
    total_reads = sum(node_reads.values())
    ax.set_title(
        f'{rep}\n{total_reads} reads,  {sub.number_of_nodes()}n / '
        + f'{sub.number_of_edges()}e',
        fontsize=title_font,
    )
    ax.set_axis_off()

    cbar = fig.colorbar(
        nodes, ax=ax, fraction=0.04, pad=0.02, shrink=0.7
    )
    cbar.ax.tick_params(labelsize=max(6, title_font - 2))
    cbar.set_label('reads / node', fontsize=max(6, title_font - 2))

    return rep, total_reads


def plot_top_clusters(graph, out_path, n_top=20, v=0, node_label_max_nodes=500,
                      edge_label_max_edges=100):
    """Plot the top-N connex components as small networkx panels.
    Nodes are labelled with their read count; edges with (alignment_distance /
    overlap_distance). Here "top-N" means components with highest sum of node
    weight.

    Arguments:
    graph (networkx graph) - Networkx graph object.
    out_path         (str) - Path to output plot file.
    n_top            (int) - Number of (biggest) graph to plot.
    v                (int) - Level of verbosity (default: 0 = muted)
    node_label_max_nodes (int) - Skip per-node read-count labels above this
                                 cluster size (default: 500).
    edge_label_max_edges (int) - Skip per-edge distance labels above this
                                 edge count (default: 100).
    """
    if v > 0:
        t_zero = time.perf_counter()
    top = _rank_clusters(graph)[:n_top]
    if not top:
        return 0

    n_cols = 4
    n_rows = (len(top) + n_cols - 1) // n_cols
    fig, axes = plt.subplots(
        n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows)
    )
    axes = np.atleast_1d(axes).flatten()
    cmap = plt.get_cmap('viridis')

    for idx, (_, sub) in enumerate(top):
        _draw_cluster(
            sub, axes[idx], fig,
            node_size=100,
            node_label_font=6,
            edge_label_font=4,
            title_font=8,
            node_label_max_nodes=node_label_max_nodes,
            edge_label_max_edges=edge_label_max_edges,
            cmap=cmap,
        )
    for idx in range(len(top), len(axes)):
        axes[idx].set_axis_off()

    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Built cluster visualization for {len(top)} clusters in '
            + f'{elapsed:.3f} s', v, 2, 1
        )
        utils.send_text('Saving clusters visualization plot', v, 3, 1)

    plt.tight_layout()
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close(fig)

    return 0


def plot_clusters_individually(graph, out_dir, n_top=20, v=0,
                               node_label_max_nodes=500,
                               edge_label_max_edges=100):
    """One PNG per top-N cluster, in their own directory. Bigger figure,
    bigger fonts, and higher label thresholds than the combined view —
    designed for dense clusters that get crushed in the panel layout.
    Files are named `cluster_<rank>_<representative>.png`.

    Arguments:
    graph (networkx graph) - Networkx graph object.
    out_dir          (str) - Path to output directory.
    n_top            (int) - Number of (biggest) graph to plot.
    v                (int) - Level of verbosity (default: 0 = muted)
    node_label_max_nodes (int) - Skip per-node read-count labels above this
                                 cluster size (default: 500).
    edge_label_max_edges (int) - Skip per-edge distance labels above this
                                 edge count (default: 100).
    """
    if v > 0:
        t_zero = time.perf_counter()
    os.makedirs(out_dir, exist_ok=True)
    top = _rank_clusters(graph)[:n_top]
    if not top:
        return 0

    cmap = plt.get_cmap('viridis')
    for rank, (_, sub) in enumerate(top, start=1):
        fig, ax = plt.subplots(figsize=(14, 11))
        rep, _ = _draw_cluster(
            sub, ax, fig,
            node_size=200,
            node_label_font=8,
            edge_label_font=6,
            title_font=12,
            node_label_max_nodes=node_label_max_nodes,
            edge_label_max_edges=edge_label_max_edges,
            cmap=cmap,
        )
        # Filesystem-safe representative name
        safe_rep = rep.replace('/', '_').replace(' ', '_')
        out_path = os.path.join(
            out_dir, f'cluster_{rank:03d}_{safe_rep}.png'
        )
        plt.tight_layout()
        plt.savefig(out_path, dpi=200, bbox_inches='tight')
        plt.close(fig)

    if v > 0:
        elapsed = time.perf_counter() - t_zero
        utils.send_text(
            f'Plotted {len(top)} individual cluster PNGs in '
            + f'{elapsed:.3f} s', v, 2, 1
        )

    return 0


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-g', '--graph_pickle', type=str, required=True,
                        help='Path to input pickle file containing the graph.')
    parser.add_argument('-o', '--output', type=str, required=True,
                        help='Path to the output PNG file containing the '
                        "graph's biggest connected (connex) components "
                        'visualization.')
    parser.add_argument('-n', '--nb_biggest', type=int, default=20,
                        help='Number of largest connected components to plot '
                        "(in terms of sum of nodes weights in the component). "
                        'Default: 20.')
    parser.add_argument('-d', '--individual_dir', type=str, default=None,
                        help='Optional directory where to write one PNG per '
                        'cluster (top-N, same ranking as the combined plot). '
                        'Bigger frames and larger fonts than the combined PNG;'
                        ' useful for dense clusters. Default: skip.')
    parser.add_argument('-v', '--verbose', type=int, default=0,
                        help='Level of verbosity (default: 0 = muted)')

    args = parser.parse_args()
    v = args.verbose

    # Ensure output directory exists
    output_dir = os.path.split(args.output)[0]
    if len(output_dir) > 0 and not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # Load graph
    utils.send_text(f'Loading graphs', v, 1, 0)
    with open(args.graph_pickle, 'rb') as graph_pkl_file:
        graph = pickle.load(graph_pkl_file)

    # Plot connected components
    utils.send_text('Plotting biggest connected components', v, 1, 0)
    res = plot_top_clusters(
        graph,
        args.output,
        n_top=args.nb_biggest,
        v=v
    )

    # Plot individual png graphs
    if args.individual_dir is not None:
        utils.send_text('Plotting individual cluster PNGs', v, 1, 0)
        plot_clusters_individually(
            graph,
            args.individual_dir,
            n_top=args.nb_biggest,
            v=v,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
