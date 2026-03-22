from typing import Union




class PyGraphVisual:
    """
    Facilitate data preparation for graph visualization that uses the Cytoscape.js library

        NOTE: this is a fork of a class by the same name in the sister project BrainAnnex.org (where the leading
          development for this module takes place)
    """


    def __init__(self, db=None):
        self.db = db                    # Object of "GraphAccess" class

        self.nodes = []                 # The node data needed by the Cytoscape graph visualization.
                                        #   A list of dicts defining nodes.
                                        #   NODES must contain the keys 'id' and '_node_labels'
                                        #   EXAMPLE (2 nodes):
                                        #   [{'id': 1, '_node_labels': ['PERSON'], 'name': 'Julian'},
                                        #    {'id': 2, '_node_labels': ['CAR'], 'color': 'white'}]

        self.edges = []                 # The edge data needed by the Cytoscape graph visualization.
                                        #   A list of dicts defining edges.
                                        #   EDGES must contain the keys 'name', 'source', and 'target'
                                        #       (and presumably 'id' is required as well?)
                                        # EXAMPLE (1 edge):
                                        #   [{'name': 'OWNS', 'source': 1, 'target': 2, 'id': 'edge-1'}]
                                        #   The above assumes that nodes 1 and 2 exist

        self.color_mapping = {}         # Mapping a node label to its interior color (the edge color is an automatic variation)
                                        # EXAMPLE:  {'PERSON': 'cyan', 'CAR': 'orange'}

        self.caption_mapping = {}       # Mapping a node label to its caption (field name to use on the graph)
                                        # EXAMPLE:  {'PERSON': 'name', 'CAR': 'color'}

        self._all_node_ids = []         # List of all the node id's added to the graph so far;
                                        #   used for optional prevention of duplicates

        self._next_available_edge_id = 1 # An auto-increment value


        self.EXTRA_COLORS =  {          # Convenient extra colors, not available thru standard CSS names
            "graph_green": '#8DCC92',
            "graph_darkgreen": '#56947E',
            "graph_teal": '#569481',
            "graph_orange": '#F79768',
            "graph_blue": '#4C8EDB',
            "graph_red": '#F16668',
            "graph_darkbrown": '#604A0D',
            "graph_lightbrown": '#D9C8AD',
            "graph_orchid": '#C990C1',
            "graph_gold": '#FFC455'     # Has a bit more orange than 'gold'
        }



    def __str__(self) -> str:
        """
        Return an overview description of this object

        :return:    A string with a description of this object
        """
        s = f"Object describing a graph containing {len(self.nodes)} nodes " \
            f"and {len(self.edges)} edges:\n"

        s += f"    Graph Node structure: {self.nodes}\n"
        s += f"    Graph Edge structure: {self.edges}\n"
        s += f"    Graph color mapping: {self.color_mapping}\n"
        s += f"    Graph caption mapping: {self.caption_mapping}"

        return s



    def get_graph_data(self) -> dict:
        """
        Return a data dictionary with all the relevant visualization info
        (typically, to pass to the front-end, for eventual use by the Cytoscape.js library)

        :return:    A dict with 4 keys - "nodes", "edges", "color_mapping" and "caption_mapping"
                        For more details, see object constructor
        """
        #TODO: perhaps add a 2nd arg, to specify just one element to be returned
        return {"nodes": self.nodes,
                "edges": self.edges,
                "color_mapping": self.color_mapping,
                "caption_mapping": self.caption_mapping}



    def add_node(self, node_id :int|str, labels="", properties=None) -> None:
        """
        From the given data about a database node, assemble and cumulatively store the data
        in a format expected by Cytoscape.js (the visualization front-end).

        Calls with a duplicate node_id are silently disregarded.

        NO database operation is performed.  You need to first extract, or have on hand, that node data.

        EXAMPLE:    given a call to:  add_node(node_id=1, labels=['PERSON'], data={'name': 'Julian'})
                    then the internal format, cumulatively added to self.nodes, will be:
                    {'id': 1, '_node_labels': ['PERSON'], 'name': 'Julian'}

        :param node_id:     Either an integer or a string to uniquely identify this node in the graph;
                                it will be used to specify the start/end nodes of edges.
                                Typically, use either the internal database ID, or some URI.
                                Records with a duplicate node_id are silently disregarded
        :param labels:      A string, or list of strings, with the node's label(s) used in the graph database;
                                if not used, pass an empty string
        :param properties:  A dict with all the node properties of interest
        :return:            None
        """
        assert isinstance(node_id, (int, str)), \
            "add_node(): the argument `node_id` must be an integer or a string"

        if type(node_id) == str:
            assert node_id != "", "add_node(): cannot use an empty string for the argument `node_id`"


        if node_id in self._all_node_ids:
            return      # Silently disregard duplicates.   TODO: use a "sorted list" for `self._all_node_ids`

        if properties is None:
            properties = {}

        d = properties.copy()     # Make an independent clone of the data dict


        d["id"] = node_id

        if type(labels) == str:
            labels = [labels]   # Turn into a list, if passed as a single string value

        d["_node_labels"] = labels

        self.nodes.append(d)
        self._all_node_ids.append(node_id)



    def add_edge(self, from_node :str|int, to_node :str|int,
                 name :str, edge_id=None, properties=None) -> None:
        """
        Prepare and store the data for 1 edge, in a format expected by the visualization front-end.

        EXAMPLE:    given a call to:  add_edge(from_node=1, to_node=2, name='OWNS', edge_id='edge-1')
                    then the internal format, cumulatively added to self.edges, will be:
                    {'name': 'OWNS', 'source': 1, 'target': 2, 'id': 'edge-1'}

        Note that the id's, in this context, is whatever we pass to the visualization module for the nodes;
        not necessarily related to the internal database ID's

        :param from_node:   Integer or a string uniquely identify the "id" of the node where the edge originates
        :param to_node:     Integer or a string uniquely identify the "id" of the node where the edge ends
        :param name:        Name of the relationship
        :param edge_id:     [OPTIONAL]  If not provided, strings such as "edge-123" are used,
                                with auto-generated consecutive integers
        :param properties:  [OPTIONAL]  A dict with all the edge properties of interest.
                                If keys with conflicting names are present ('name', 'source', 'target', 'id'),
                                an underscore is prefixed to them

        :return:            None
        """
        # TODO: check if the edge already exists, with same name
        # TODO: unclear if edge id's are really needed by Cytoscape.js

        d = {"name": name, "source": from_node, "target": to_node}

        if edge_id is not None:
            d["id"] = edge_id
        else:
            # Use an auto-increment value of the form "edge-n" for some n
            d["id"] = f"edge-{self._next_available_edge_id}"
            self._next_available_edge_id += 1        # Maintain an auto-increment value for edge ID's

        if properties:
            if "name" in properties:
                properties["_name"] = properties["name"]
                del properties["name"]
            if "source" in properties:
                properties["_source"] = properties["source"]
                del properties["source"]
            if "target" in properties:
                properties["_target"] = properties["target"]
                del properties["target"]
            if "id" in properties:
                properties["_id"] = properties["id"]
                del properties["id"]

            d.update(properties)    # Also store the edge properties in its dict

        self.edges.append(d)



    def assign_caption(self, label :str, caption="name") -> None:
        """
        Incrementally assign and store a mapping from label name to caption (name of field to use on the display).
        Additional calls to this function will register further assignments.

        EXAMPLES:   assign_caption(label='PERSON', caption='name')
                    assign_caption(label='CAR', caption='color')

        :param label:   The name of a node label
        :param caption: The name of the node field to use on the display (by default, "name")
        :return:        None
        """
        assert type(label) == str, \
            "assign_caption(): the argument `label` must be a string"

        self.caption_mapping[label] = caption



    def assign_color_mapping(self, label :str, color :str) -> None:
        """
        Incrementally assign and store a mapping from label name to color to use for nodes having that label.
        Additional calls to this function will register further assignments.

        EXAMPLES:   assign_color_mapping(label='PERSON', color='purple')
                    assign_color_mapping(label='CAR', color='#FF0088')

        :param label:   The name of a node label
        :param color:   The name (hex string or standard CSS name) of the color to use on the display,
                            OR one of the non-standard names provided by this class,
                            in the object variable self.EXTRA_COLORS
        :return:        None
        """
        if color in self.EXTRA_COLORS:
            color = self.EXTRA_COLORS[color]     # Convert the non-standard name into a numerical value

        self.color_mapping[label] = color




    ############   The methods below require a database connection   ############


    def prepare_recordset(self, id_list :[int|str]) -> [dict]:
        """
        Given a list of node internal id's, construct and return a dataset of their properties,
        plus the special fields `internal_id` and `_node_labels`

        :param id_list: A list of internal id's of database nodes

        :return:        A list of dicts, with the node properties plus the special fields `internal_id` and `_node_labels`
                            EXAMPLE:
                            [{'_internal_id': 123, '_node_labels': ['Person'], 'name': 'Julian'}]
        """
        # TODO: move to GraphAccess
        assert type(id_list) == list, \
            f"prepare_recordset(): argument `id_list` must be a list; it is of type {type(id_list)}"

        assert self.db, \
            "prepare_recordset(): missing database handle; did you pass it when instantiating PyGraphVisual(db=...) ?"

        if id_list == []:
            return []       # No data was passed

        q = f'''
            MATCH (n)
            WHERE ID(n) IN $id_list
            RETURN n
            '''

        #self.db.debug_print_query(q, {"id_list": id_list})
        return self.db.query_extended(q, {"id_list": id_list}, flatten=True)



    def prepare_graph(self, result_dataset :[dict], cumulative=False, add_edges=True, avoid_links=None) -> [int|str]:
        """
        Given a list of dictionary data containing the properties of graph-database nodes - for example,
        as returned by GraphAccess.get_nodes() - construct, and save inside this pyhon object,
        the visualization data for them.

        Each passed dictionary entry MUST have a key named "_internal_id".
        If any key named "id" is found, it get automatically renamed "_id_original" (since "id" is used by the visualization software);
        if "_id_original" already exists, an Exception is raised.  (Copies are made; the original data object isn't affected.)
        Though not required, a key named "_node_labels" is typically present as well.

        Any date/datetime value found in the database will first be "sanitized" into a string representation of the date;
        the time portion, if present, will get dropped

        :param result_dataset:  A list of dictionary data about graph-database nodes;
                                    each dict must contain an entry with the key "_internal_id".
                                    No problem if duplicates in "_internal_id" are present;
                                    only the fist of one duplicates gets processed
        :param cumulative:      If False (default) then any previous call to this function will get ignored,
                                    and a new graph is appended
        :param add_edges:       If True, all existing edges among the displayed nodes
                                    will also be part of the visualization
        :param avoid_links:     Name or list of name of links to avoid including

        :return:                A list of values with the internal databased IDs
                                    of all the nodes added to the graph structure
        """
        assert self.db, \
            "prepare_graph(): missing database handle; did you pass it when instantiating PyGraphVisual(db=...) ?"

        assert type(result_dataset) == list, \
            f"prepare_graph(): argument `id_list` must be a list; it is of type {type(result_dataset)}"

        if len(result_dataset) == 0:
            return []       # No data was passed


        if not cumulative:
            # Reset: clear out any previous graph-structure data, and reset edge number auto-increment
            self.nodes = []
            self.edges = []
            self._all_node_ids = []
            self._next_available_edge_id = 1


        id_key_renaming = False
        node_list = []      # Running list of internal databased IDs, for nodes in `result_dataset`
        for node in result_dataset:
            internal_id = node.get("_internal_id")
            assert internal_id is not None, \
                f"prepare_graph() - the following record lacks the required `internal_id` key: {node}"

            if internal_id in node_list:
                # Ignore records already seen (duplicates in `internal_id`)
                print(f"prepare_graph(): ignoring duplicate record with `internal_id` = {internal_id}")
                continue

            node_clone = node.copy()

            if "id" in node_clone:
                assert "_id_original" not in node_clone, \
                    f"prepare_graph(): keys named `id` are routinely automatically renamed `_id_original`, " \
                    f"but the latter key also already exists!  Found in: {node_clone}"
                node_clone["_id_original"] = node_clone["id"]
                del node_clone["id"]
                id_key_renaming = True

            node_list.append(internal_id)

            if "_node_labels" in node_clone:
                labels = node["_node_labels"]
            else:
                labels = ""

            self.add_node(node_id=internal_id, labels=labels,
                          properties=self.db.sanitize_date_times(node_clone, drop_time=True))


        if add_edges:
            # Search the database for any edges among any of the nodes selected for the visualization
            exclude_clause = ""
            if avoid_links:
                if type(avoid_links) == str:
                    avoid_links = [avoid_links]
                else:
                    assert type(avoid_links) == list,  \
                        f"prepare_graph(): The argument `avoid_links`, if passed, must be a string or a list"

                exclude_clause = f"AND NOT type(r) IN {avoid_links}"


            q = f'''
                MATCH (n1)-[r]->(n2) 
                WHERE ID(n1) IN $node_list AND ID(n2) IN $node_list 
                {exclude_clause}
                RETURN DISTINCT id(n1) AS from_node, id(n2) AS to_node, 
                       type(r) AS rel_name, properties(r) AS rel_props
                '''

            #self.db.debug_print_query(q, {"node_list": node_list})
            result = self.db.query(q, {"node_list": node_list})


            # "Sanitize" the records, as needed, i.e. make them suitable for JSON serialization,
            # in anticipation of eventually passing the data to JavaScript
            for edge in result:
                #print(edge)
                edge_props = self.db.sanitize_date_times(edge["rel_props"], drop_time=True)
                self.add_edge(from_node=edge["from_node"], to_node=edge["to_node"],
                              name=edge["rel_name"], properties=edge_props)

        if id_key_renaming:
            print("prepare_graph(): keys named `id` were found in one or more of the records; "
                  "they were renamed `_id_original` to avoid conflict with internal database IDs")

        return node_list



    def assemble_graph(self, id_list :[int|str]) -> ([dict],[dict]):
        """
        Given a list of node internal id's, construct and return
        the data needed by the Cytoscape graph visualization (including all properties).

        Any date/datetime value found in the database will first be "sanitized"
        into a string representation of the date;
        the time portion, if present, will get dropped

        :param id_list: A list of internal id's of database nodes
        :return:        A pair with all the data needed by the Cytoscape graph visualization:
                            1) list of dicts defining nodes
                            2) list of dicts defining edges
        """
        recordset = self.prepare_recordset(id_list)
        self.prepare_graph(result_dataset=recordset, cumulative=False, add_edges=True)
        return (self.nodes , self.edges)



    def link_node_groups(self, group1 :[int|str], group2 :[int|str]) -> None:
        """
        Search the database for any edges from any of the nodes in the 1st group, to any node in the 2nd group.
        Any located edge will be added to the visualization data stored in this object.

        :param group1:  List of the internal databased IDs of the 1st group of database nodes
        :param group2:  List of the internal databased IDs of the 2nd group of database nodes
        :return:        None
        """
        q = '''
            MATCH (n1)-[r]->(n2) 
            WHERE ID(n1) IN $group1 AND ID(n2) IN $group2 
            RETURN DISTINCT id(n1) AS from_node, id(n2) AS to_node, type(r) AS name
            '''

        result = self.db.query(q, {"group1": group1, "group2": group2})
        for edge in result:
            #print(edge)
            self.add_edge(from_node=edge["from_node"], to_node=edge["to_node"], name=edge["name"])



##########################################################################################

class PyGraphVisual_OLD:
    """
    Facilitate data preparation for graph visualization using the Cytoscape.js library.

    NOTE: this is a fork of a class by the same name in the sister project BrainAnnex.org (where the leading
          development for this module takes place)
    """


    def __init__(self, db=None):
        self.db = db                    # Object of "GraphAccess" class

        self.structure = []             # The data that defines a graph to visualize.
                                        #   A list of dicts defining nodes, and possibly edges as well.
                                        #   NODES must contain the keys 'id' and 'labels'
                                        #   EDGES must contain the keys 'name', 'source', and 'target'
                                        #       (and presumably 'id' is required as well?)
                                        # EXAMPLE (2 nodes and an edge):
                                        #   [{'id': 1, 'labels': ['PERSON'], 'name': 'Julian'},
                                        #    {'id': 2, 'labels': ['CAR'], 'color': 'white'},
                                        #    {'name': 'OWNS', 'source': 1, 'target': 2, 'id': 'edge-1'}]

        self.color_mapping = {}         # Mapping a node label to its interior color (the edge color is an automatic variation)
                                        # EXAMPLE:  {'PERSON': 'cyan', 'CAR': 'orange'}

        self.caption_mapping = {}       # Mapping a node label to its caption (field name to use on the graph)
                                        # EXAMPLE:  {'PERSON': 'name', 'CAR': 'color'}

        self._next_available_edge_id = 1 # An auto-increment value

        self._all_node_ids = []         # List of all the node id's added to the graph so far;
                                        #   used for optional prevention of duplicates


    def serialize(self) -> dict:
        """
        Extract and return a dict of relevant data from this object

        :return:
        """
        return {"structure": self.structure,
                "color_mapping": self.color_mapping,
                "caption_mapping": self.caption_mapping}



    def __str__(self) -> str:
        """
        Return an overview description of this object

        :return:    A string with a description of this object
        """
        s = f"Object describing a graph containing a total of {len(self.structure)} nodes and edges:\n"
        s += f"    Graph structure: {self.structure}\n"
        s += f"    Graph color mapping: {self.color_mapping}\n"
        s += f"    Graph caption mapping: {self.caption_mapping}"

        return s



    def get_graph_data(self) -> dict:
        """
        Return a data dictionary with all the relevant visualization info
        (typically, to pass to the front-end, for eventual use by the Cytoscape.js library)

        :return:    A dict with 3 keys - "structure", "color_mapping" and "caption_mapping"
        """
        return {"structure": self.structure,
                "color_mapping": self.color_mapping,
                "caption_mapping": self.caption_mapping}



    def add_node(self, node_id :Union[int,str], labels="", data=None) -> None:
        """
        Prepare and store the data for 1 node, in a format expected by the visualization front-end

        EXAMPLE:    {'id': 1, 'name': 'Julian', 'labels': ['PERSON']}

        :param node_id:         Either an integer or a string to uniquely identify this node in the graph;
                                    it will be used to specify the start/end nodes of edges.
                                    Typically, use either the internal database ID, or some URI
                                    Records with a duplicate node_id are silently disregarded
        :param labels:          A string, or list of strings, with the node's label(s) used in the graph database;
                                    if not used, pass an empty string
        :param data:            A dict with all other node data not already specified in any of the other arguments
        :return:                None
        """
        assert node_id != "", \
            "add_node(): cannot use an empty string for the argument `node_id`"


        if node_id in self._all_node_ids:
            return      # Silently disregard duplicates

        if data is None:
            data = {}

        d = data.copy()     # Make an independent clone of the data dict


        d["id"] = node_id

        if type(labels) == str:
            labels = [labels]

        d["labels"] = labels

        self.structure.append(d)
        self._all_node_ids.append(node_id)



    def add_edge(self, from_node :Union[str, int], to_node :Union[str, int],
                 name :str, edge_id=None, data=None) -> None:
        """
        Prepare and store the data for 1 edge, in a format expected by the visualization front-end.

        EXAMPLE:   {'name': 'OWNS', 'source': 1, 'target': 2, 'id': 'edge-1'}

        Note that 'id', in this context, is whatever we pass to the visualization module for the nodes;
        not necessarily related to the internal database ID's

        :param from_node:   Integer or a string uniquely identify the "id" of the node where the edge originates
        :param to_node:     Integer or a string uniquely identify the "id" of the node where the edge ends
        :param name:        Name of the relationship
        :param edge_id:     (OPTIONAL)  If not provided, strings such as "edge-123" are used,
                                        with auto-generated consecutive integers
                                TODO: unclear if edge id's are really needed
        :param data:        A dict with all other node data not already specified in any of the other arguments

        :return:            None
        """
        from_id = from_node
        to_id = to_node

        d = {"name": name, "source": from_id, "target": to_id}

        if edge_id is not None:
            d["id"] = edge_id
        else:
            d["id"] = f"edge-{self._next_available_edge_id}"
            self._next_available_edge_id += 1        # Maintain an auto-increment value

        if data:
            d.update(data)

        self.structure.append(d)



    def assign_caption(self, label :str, caption="name") -> None:
        """
        Assign and store a mapping from label name to caption (name of field to use on the display)

        EXAMPLES:   assign_caption(label='PERSON', caption='name')
                    assign_caption(label='CAR', caption='color')

        :param label:   The name of a node label
        :param caption: The name of the node field to use on the display (by default, "name")
        :return:        None
        """
        assert type(label) == str, \
            "assign_caption(): the argument `label` must be a string"

        self.caption_mapping[label] = caption



    def assign_color_mapping(self, label :str, color :str) -> None:
        """
        Assign and store a mapping from label name to color to use for nodes having that label

        :param label:   The name of a node label
        :param color:   The name (hex string or standard CSS name) of the color to use on the display
        :return:        None
        """
        extra_colors =  {       # Convenient extra colors, not available thru standard CSS names
            "graph_green": '#8DCC92',
            "graph_darkgreen": '#56947E',
            "graph_teal": '#569481',
            "graph_orange": '#F79768',
            "graph_blue": '#4C8EDB',
            "graph_red": '#F16668',
            "graph_darkbrown": '#604A0D',
            "graph_lightbrown": '#D9C8AD',
            "graph_orchid": '#C990C1',
            "graph_gold": '#FFC455'     # Has a bit more orange than 'gold'
        }

        if color in extra_colors:
            color = extra_colors[color]

        self.color_mapping[label] = color




    ##########   The method below require a database connection   ##########


    def prepare_graph(self, result_dataset :[dict], add_edges=True) -> [int]:
        """
        Given a list of dictionary data about graph-database nodes - for example,
        as returned by GraphAccess.get_nodes() - construct and save visualization data for them.

        Each dictionary entry is expected to have a key named "internal_id";
        if not present, it will be silently ignored.
        Though not required, a key named "neo4j_labels" is typically present as well.

        :param result_dataset:  A list of dictionary data about graph-database nodes
        :param add_edges:       If True, all existing edges among the displayed nodes
                                    will also be part of the visualization
        :return:                A list of integers with the internal databased IDs
                                    of all the nodes added to the graph structure
        """
        node_list = []      # Running list of internal databased IDs, for nodes included in visualization
        for node in result_dataset:
            internal_id = node.get("internal_id")
            if not internal_id:
                continue    # Silently ignore any entry lacking this data

            node_list.append(internal_id)

            if "neo4j_labels" in node:
                labels = node["neo4j_labels"]
                del node["neo4j_labels"]
            else:
                labels = ""

            self.add_node(node_id=internal_id, labels=labels,
                          data=node)


        if add_edges:
            # Search the database for any edges among any of the nodes selected for the visualization
            q = '''
                MATCH (n1)-[r]->(n2) 
                WHERE ID(n1) IN $node_list AND ID(n2) IN $node_list 
                RETURN DISTINCT id(n1) AS from_node, id(n2) AS to_node, type(r) AS name
                '''

            result = self.db.query(q, {"node_list": node_list})
            for edge in result:
                #print(edge)
                self.add_edge(from_node=edge["from_node"], to_node=edge["to_node"], name=edge["name"])

        return node_list



    def link_node_groups(self, group1 :[int], group2 :[int]) -> None:
        """
        Search the database for any edges from any of the nodes in the 1st group, to any node in the 2nd group.
        Any located edge will be added to the visualization data stored in this object

        :param group1:  List of integers with the internal databased IDs of the 1st group of nodes
        :param group2:  List of integers with the internal databased IDs of the 2nd group of nodes
        :return:        None
        """
        q = '''
            MATCH (n1)-[r]->(n2) 
            WHERE ID(n1) IN $group1 AND ID(n2) IN $group2 
            RETURN DISTINCT id(n1) AS from_node, id(n2) AS to_node, type(r) AS name
            '''

        result = self.db.query(q, {"group1": group1, "group2": group2})
        for edge in result:
            #print(edge)
            self.add_edge(from_node=edge["from_node"], to_node=edge["to_node"], name=edge["name"])
