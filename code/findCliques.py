import argparse
import pandas as pd
import networkx as nx

def filter_fully_connected_components(input_file, output_node_file, output_edge_file, min_clique_size):
    # Read the input CSV file
    df = pd.read_csv(input_file)
    df = df.drop(columns=['Annotation', 'Score', 'EdgeType']).drop_duplicates()

    print("Initial number of rows:", len(df))
    
    # Create the graph
    G = nx.from_pandas_edgelist(df, 'ID1', 'ID2')

    # Find all maximal cliques in the graph
    cliques = list(nx.find_cliques(G))
    
    # Initialize lists for node table and edge table
    node_table = []
    edge_table = []

    # Create a dictionary to store the largest clique ID for each node
    node_clique_status = {node: None for node in G.nodes()}
    # Create a set to store edges that are part of cliques
    edges_in_cliques = set()
    # Create a dictionary to store the final cliques
    final_cliques = {}

    # Assign clique IDs and mark nodes and edges
    for clique_id, clique in enumerate(cliques):
        clique_size = len(clique)
        formatted_clique_id = f"{clique_id}_{clique_size}"
        
        # Handle cliques with two members separately
        if clique_size == 2:
            node1, node2 = clique
            if len(list(G.neighbors(node1))) == 1 and len(list(G.neighbors(node2))) == 1:
                final_cliques[formatted_clique_id] = set(clique)
                for node in clique:
                    current_clique_size = int(node_clique_status[node].split('_')[1]) if node_clique_status[node] else 0
                    if clique_size > current_clique_size:
                        node_clique_status[node] = formatted_clique_id
                edges_in_cliques.add(tuple(sorted((node1, node2))))
        elif clique_size >= min_clique_size:
            final_cliques[formatted_clique_id] = set(clique)
            for node in clique:
                current_clique_size = int(node_clique_status[node].split('_')[1]) if node_clique_status[node] else 0
                if clique_size > current_clique_size:
                    node_clique_status[node] = formatted_clique_id
            for i in range(clique_size):
                for j in range(i + 1, clique_size):
                    edge = tuple(sorted((clique[i], clique[j])))
                    edges_in_cliques.add(edge)

    # Reassign nodes to their final cliques
    for node, clique_id in node_clique_status.items():
        if clique_id:
            final_cliques[clique_id].add(node)

    # Update clique IDs to reflect final sizes
    new_clique_ids = {}
    for old_id, nodes in final_cliques.items():
        new_id = f"{old_id.split('_')[0]}_{len(nodes)}"
        for node in nodes:
            node_clique_status[node] = new_id
        new_clique_ids[old_id] = new_id

    # Filter edges to only keep those within final cliques
    final_edges = set()
    for edge in edges_in_cliques:
        node1_clique = node_clique_status[edge[0]]
        node2_clique = node_clique_status[edge[1]]
        if node1_clique == node2_clique:
            final_edges.add(edge)

    # Save node attributes to the node table
    nodes_in_cliques = 0
    for node, clique_id in node_clique_status.items():
        if clique_id:  # Only include nodes that are part of at least one clique
            node_table.append({'Node': node, 'clique_id': clique_id})
            nodes_in_cliques += 1

    # Save edge attributes to the edge table
    for edge in final_edges:
        edge_table.append({'Source': edge[0], 'Target': edge[1], 'clique': True})

    # Convert node and edge tables to dataframes
    node_df = pd.DataFrame(node_table)
    edge_df = pd.DataFrame(edge_table)
    
    # Write the node and edge tables to CSV files
    node_df.to_csv(output_node_file, index=False)
    edge_df.to_csv(output_edge_file, index=False)

    # Print statistics
    total_nodes = len(G.nodes())
    nodes_not_in_cliques = total_nodes - nodes_in_cliques
    print(f"Total cliques found: {len(final_cliques)}")
    print(f"Total nodes in cliques: {nodes_in_cliques}")
    print(f"Total nodes not in cliques: {nodes_not_in_cliques}")

    print("Final number of node rows:", len(node_df))
    print("Final number of edge rows:", len(edge_df))

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Filter fully connected components and add clique attributes.')
    parser.add_argument('input_file', type=str, help='Path to the input CSV file.')
    parser.add_argument('output_node_file', type=str, help='Path to the output node CSV file.')
    parser.add_argument('output_edge_file', type=str, help='Path to the output edge CSV file.')
    parser.add_argument('--min_clique_size', type=int, default=2, help='Minimum size of cliques to be considered.')
    
    args = parser.parse_args()
    
    filter_fully_connected_components(args.input_file, args.output_node_file, args.output_edge_file, args.min_clique_size)
