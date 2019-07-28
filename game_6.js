// game principle
/*
circular maze, a central circle (player) needs to get out in time
problem: circle grows - when too large, player dies (crushed between maze walls)
circle size can be resetted when collecting items ("meds")
player control via arrow keys
*/

// init global vars
var n_rings;
var ring_spacing;
var ring_thickness;
var n_meds;
var size_grow;
var player_size;
var player_speed;
var collision_steps;
var max_holes;
var dead_end_frac;
var max_distance;
var med_min_distance;
var med_max_distance;
var rings = [];
var walls = [];
var start_rot;
var hole_problem_level = 1; // layer at which overlap problems need to be resolved

// debug flags
var no_virt_holes = false;
var debug_draw_toggle = false;

// add event listenersggg
document.addEventListener("keydown", keydown);
document.addEventListener("keyup", keyup);
document.addEventListener("click", mouseclick);

// prepare values for vars
function add_values() {
    // global vars that are initialized before
    n_rings = 3;            // rings of maze
    ring_spacing = 80;      // spacing of maze rings in px
    ring_thickness = 15;    // diameter of maze walls in px
    n_meds = 5;             // how many items need to be collected before allowed to exit
    size_grow = 0.005;          // increase of player diameter over a frame
    player_size = 3;        // initial player diameter in px
    player_speed = 5;       // movement speed
    collision_steps = 10;   // determines how many times per draw the collision resolution function is called
    max_holes = 6;          // number of holes for middle rings
    dead_end_frac = 0.2;       // percentage of nodes in typical layer that are dead ends
    med_min_distance_fac = 0.2;  // factor of max distance that player is separated by med
    med_max_distance_fac = 0.8;  // how far is last med from player?
    start_rot = Math.random(); // rotation of innermost ring
    // derived vars
    // TODO: max distance   // formula that takes into account player_speed, size_grow and ring_spacing
    // max_distance = ;
}

// classes

class MazeNetwork {

    // TODO layers are lists of nodes, i.e. network --> replace layers with network
    constructor(n_layers=n_rings, max_open_nodes=max_holes) {
        // input class variables
        this.n_layers = n_layers;
        this.max_open_nodes = max_open_nodes;

        // setup network (random node positions)
        this.n_nodes_network = this.generate_number_holes();
        this.n_dead_ends_network = this.generate_dead_end_numbers();
        this.n_virtual_nodes_network = this.generate_virtual_numbers();

        // create first network layout and loop randomly until valid
        // valid means: no overlaps in node order
        this.create_network_layout();
        while (this.overlap_network()) {
            this.create_network_layout();
        }

    }
    create_network_layout() {
        this.layers = this.setup_layers_with_nodes();
        this.setup_dead_ends();
        this.add_virt_nodes();
        this.connect_nodes();
    }
    generate_number_holes() {

        // following schmematic: 1,2,2,4,4,8,8,16,16,...,1
        var num = 2;
        var num_seq = [1];
        while (num_seq.length < this.n_layers - 1) {
            num_seq.push(num);
            num_seq.push(num);
            if (num < 8) { num *= 2; } // clip to upper bound
        }
        // handle case if even number --> add last num
        if (num_seq.length < this.n_layers) {
            num_seq.push(num);
        }
        num_seq[num_seq.length-1] = 1;
        return num_seq;
    }
    generate_dead_end_numbers() {
    
        // 0,1,0,1,0,2,0,4,0,8,0,...
        var num = 1;
        var num_seq = [];
        while (num_seq.length < this.n_layers - 1) {
            num_seq.push(0);
            num_seq.push(num);
            if (num_seq.length >= 4) {
                // clip to upper bound
                if (num < 2) {
                    num *= 2;
                } else {
                    num *= 0;
                }
            }
        }
        // handle case if even number --> add last num
        if (num_seq.length < this.n_layers) {
            num_seq.push(0);
        }
        // last layer: no dead ends
        num_seq[num_seq.length-1] = 0;
        
        return num_seq;
    }
    generate_virtual_numbers() {
        var num = 1;
        // assuming n >= 2
        var num_seq = [0,0];
        while (num_seq.length < this.n_layers - 1) {
            num_seq.push(num);
            num_seq.push(num);
            // clip to upper bound
            if (num < 4) {
                num *= 2;
            }
        }
        // handle case if even number --> add last num
        if (num_seq.length < this.n_layers) {
            num_seq.push(num);
        }
        return num_seq;
    }
    setup_layers_with_nodes() {
        var layers = [];
        for (let index = 0; index < this.n_layers; index++) {
            layers.push([]);
            if (index == 0 || index == this.n_layers - 1) {
                // only one node
                layers[layers.length-1].push(new Node(index, 0, this));
            } else {
                var n_nodes = this.n_nodes_network[index] - this.n_virtual_nodes_network[index];
                for (let index2 = 0; index2 < n_nodes; index2++) {
                    layers[layers.length-1].push(new Node(index, index2, this));
                }
            }
        }
        return layers;
    }
    setup_dead_ends() {
        // set the dead ends, starting from last layer moving upwards
        var last_layer_size = this.max_open_nodes; // keeping track of n nodes in last layer
        for (let index = this.layers.length-1; index >= 0; index--) {
            const layer = this.layers[index];
            var n_dead_ends = this.n_dead_ends_network[index];
            // adjust if last_layer_size smaller than current layer.length
            if (last_layer_size < layer.length && n_dead_ends < last_layer_size) {
                n_dead_ends = layer.length - last_layer_size;
            }
            // different case second-last layer: all but one are dead ends
            if (index == this.layers.length-2) {
                n_dead_ends = layer.length - 1;
            }
            var shuffle_array = [];
            for (let i = 0; i < n_dead_ends; i++) {
                shuffle_array.push(true);
            }
            for (let i = 0; i < layer.length - n_dead_ends; i++) {
                shuffle_array.push(false);
            }
            shuffle_array = shuffle(shuffle_array);
            for (let i = 0; i < shuffle_array.length; i++) {
                this.layers[index][i].dead_end = shuffle_array[i];
            }
            last_layer_size = layer.length;
        }
    }
    add_virt_nodes_layer(index, index2) {
        // add virtual node in in next layer if dead end or already virtual
        if (this.layers[index][index2].dead_end || this.layers[index][index2].virtual) {
            var virt_node = new Node(index + 1, index2, this);
            virt_node.virtual = true;
            this.layers[index + 1].splice(index2, 0, virt_node);
            // refresh all following IDs in this layer
            for (let index3 = index2 + 1; index3 < this.layers[index + 1].length; index3++) {
                this.layers[index + 1][index3].id++;
            }
            this.layers[index + 1][index2].set_parent(this.layers[index][index2]);
        }
    }
    add_virt_nodes() {
        for (let index = 0; index < this.layers.length - 1; index++) {
            for (let index2 = 0; index2 < this.layers[index].length; index2++) {
                this.add_virt_nodes_layer(index, index2);
            }
        }
    }
    connect_nodes() {

        for (let index = this.layers.length-1; index > 0; index--) {
            // going backwards, stopping at second layer
            // setting the parent of each node
            var upper_passage_count = 0;
            for (let i = 0; i < this.layers[index - 1].length; i++) {
                const node = this.layers[index - 1][i];
                if (!node.dead_end && !node.virtual) { upper_passage_count++; }
            }
            // identify virtual nodes in this layer
            var nonvirtual_count = this.layers[index].length;
            var nonvirtual_ids = [];
            for (let i = 0; i < this.layers[index].length; i++) {
                const node = this.layers[index][i];
                if (node.virtual) {
                    nonvirtual_count--;
                } else {
                    nonvirtual_ids.push(node.id);
                }
            }
        
            if (nonvirtual_count == upper_passage_count) { // one child per parent
                var upper_id = 0;
                for (let i = 0; i < nonvirtual_ids.length; i++) {
                    const nonvirtual_id = nonvirtual_ids[i];
                    for (let j = upper_id; j < this.layers[index - 1].length; j++) {
                        if (!this.layers[index - 1][j].dead_end && !this.layers[index - 1][j].virtual) {
                            this.layers[index][nonvirtual_id].set_parent(this.layers[index - 1][j]);
                            upper_id = j + 1;
                            break;
                        }
                    }
                }
        
            }
            else { // multiple children per parent
        
                var n_children_exceed = nonvirtual_count - upper_passage_count;
                // pick randomly from among parents which have multiple children - store number
                for (var i=0; i<n_children_exceed; i++) {
                    var rand_id = Math.floor(Math.random()*this.layers[index - 1].length);
                    // make sure no dead end or virtual node is targetted
                    var count = 0;
                    while ((this.layers[index - 1][rand_id].dead_end || this.layers[index - 1][rand_id].virtual) && count < 1000) {
                        var rand_id = Math.floor(Math.random()*this.layers[index - 1].length);
                        count++;
                    }
                    if (count >= 1000) {
                        console.log('problem here   ' + count);
                    }
                    this.layers[index - 1][rand_id].n_children++;
                }
                var lower_id = 0;
        
                for (let i = 0; i < this.layers[index - 1].length; i++) {
                    const upper_node = this.layers[index - 1][i];
                    var n_c = upper_node.n_children + 1; // + 1 because only exceeded children are noted --> one at least
                    if (!upper_node.dead_end && !upper_node.virtual) {
                        for (let j = 0; j < n_c; j++) {
                            while (this.layers[index][lower_id].virtual) {
                                lower_id++;
                            }
                            this.layers[index][lower_id].set_parent(this.layers[index - 1][i]);
                            lower_id++;
                        }
                    }
                }
        
            }
        }
    }
    overlap_network() {
        // returns false when there is no overlap
        for (let index = 0; index < this.layers.length - 1; index++) {
            for (let index2 = 0; index2 < this.layers[index].length - 1; index2++) {
                const node = this.layers[index][index2];
                // get rightmost child
                var right_id = node.get_id_of_rightmost_child();
                for (let index3 = index2 + 1; index3 < this.layers[index].length; index3++) {
                    const node_right = this.layers[index][index3];
                    // get leftmost child
                    var left_id = node_right.get_id_of_leftmost_child();
                    if (left_id < right_id) {
                        return true;
                    }
                }
            }
        }
        // no overlap problems
        return false;
    }
    render(current_node = [-1, -1]) {
        // makes sure the algo for the logical layout of maze is correct
        function get_coord(layer_ind, node_id) {
            var y_num = 10 + layer_ind*y_increase;
            var x_num = 10 + node_id*x_increase;
            return {x: x_num, y: y_num};
        }
        var y_increase = canvas.height/(this.layers.length*1.5) - 5;
        var x_increase = canvas.width/(this.max_open_nodes*2) - 5;
        var current_y = 10;
        var current_x = 10;
        var start_x = current_x;
        for (let index = 0; index < this.layers.length; index++) {
            const l = this.layers[index];
            current_x = start_x;
            for (let index2 = 0; index2 < l.length; index2++) {
                const n = l[index2];
                var n_color = "white";
                if (n.dead_end) { n_color = "red"; }
                if (n.virtual) { n_color = "grey"; }
                if (n.ind == current_node[0] && n.id == current_node[1]) { n_color = "green"};
                draw_circ_outline_debug(5, {x: current_x, y: current_y}, "black", n_color);
                // connect to parent
                for (let index3 = 0; index3 < n.parent.length; index3++) {
                    const p = n.parent[index3];
                    draw_line_debug([get_coord(p.ind, p.id), {x: current_x, y: current_y}], "black")
                }
                // set next position
                current_x += x_increase;
            }
            current_y += y_increase;
        }
    }
}

class Node {
    // virtual unit representing the maze structure
    constructor(ind, id, network) {
        this.ind = ind;
        this.id = id;
        this.parent = [];
        this.children = [];
        this.dead_end = false;
        this.virtual = false;
        this.n_children = 0;
        this.coord = 0;
        this.network = network;
    }
    set_parent(parent) {
        this.parent = [parent];
        // automatically sets children accordingly
        // but first search if child already included
        var child_included = false;
        for (let index = 0; index < this.network.layers[parent.ind][parent.id].children.length; index++) {
            var child = this.network.layers[parent.ind][parent.id].children[index];
            if (child.ind == this.ind && child.id == this.id) {
                child_included = true;
                break;
            }
        }
        if (!child_included) {
            this.network.layers[parent.ind][parent.id].children.push(this);
            // important flag to store number of children
            this.network.layers[parent.ind][parent.id].n_children = this.network.layers[parent.ind][parent.id].children.length;
        }
    }
    get_id_of_leftmost_child() {
        var min = Infinity;
        for (let index = 0; index < this.children.length; index++) {
            var child = this.children[index];
            if (child.id < min) {
                min = child.id;
            }
        }
        return min;
    }
    get_id_of_rightmost_child() {
        var max = -Infinity;
        for (let index = 0; index < this.children.length; index++) {
            var child = this.children[index];
            if (child.id > max) {
                max = child.id;
            }
        }
        return max;
    }
}

class Player {
    constructor() {
        // appearance
        this.pos = {x: canvas.width/2, y: canvas.height/2};
        this.speed = player_speed;
        this.radius = player_size;
        this.color = "green";
        // states
        this.leftmove = false;
        this.rightmove = false;
        this.upmove = false;
        this.downmove = false;
        // debug movements
        this.CWmove = false;
        this.CCWmove = false;
        // for collision
        this.vel = {x: 0, y: 0};
        this.current_ring_ind = 0; // starts in middle
        this.collision_balls = [];
    }
    set_current_ring_ind() {
        // function of distance to center, pos and ring_spacing (ring_thickness?)
        var dist = distance(this.pos, center_coord);
        // set value of current ring ind
        this.current_ring_ind = Math.floor(dist/ring_spacing);
    }
    set_relevant_coll_balls_ring(rad, ring_ind, dist) {
        var relevant_holes_coords = maze.rings[ring_ind].holes_coords;
        // TODO write function that it makes sense
        var hole_overlap = overlap_hole_region(rad, relevant_holes_coords, false, false);
        if (hole_overlap) {
            var hole_coords = overlap_hole_region(rad, relevant_holes_coords, true, false)[1];
            // console.log(hole_coords);
            var coord1 = get_exact_coord(hole_coords.start, dist, false);
            var coord2 = get_exact_coord(hole_coords.end, dist, false);
            this.collision_balls.push(coord1);
            this.collision_balls.push(coord2);
            // console.log("OVERLAP");
        } else {
            var coll_ball_pos = get_exact_coord(rad, dist, false);
            this.collision_balls.push(coll_ball_pos);
            // console.log("NO OVERLAP");
        }
    }
    set_relevant_coll_balls_wall(player_dist_center) {
        for (let index = 0; index < maze.rings[this.current_ring_ind].walls.length; index++) {
            const wall = maze.rings[this.current_ring_ind].walls[index];
            // convert to coords
            if (!wall.dug) {
                this.collision_balls.push(get_exact_coord(wall.rad, player_dist_center, false));
            }
        }
    }
    set_relevant_coll_balls() {
        // deriving current ring ind
        this.set_current_ring_ind();
        // reset coll balls
        this.collision_balls = [];
        // first set coll balls along ring
        var inner_dist = this.current_ring_ind*ring_spacing - 0.5*ring_thickness;
        var outer_dist = (this.current_ring_ind + 1)*ring_spacing - 0.5*ring_thickness;
        var player_dist_center = distance(this.pos, center_coord);
        var player_rad = coord_to_rad(this.pos);
        player_rad = flip_rad(player_rad); // TODO find out why necessary
        // ignore the innermost pseudoring (-1)
        if (this.current_ring_ind > 0 && this.current_ring_ind < n_rings) {
            this.set_relevant_coll_balls_ring(player_rad, this.current_ring_ind - 1, inner_dist);
            this.set_relevant_coll_balls_wall(player_dist_center);
        }
        if (this.current_ring_ind < n_rings) {
            this.set_relevant_coll_balls_ring(player_rad, this.current_ring_ind, outer_dist);
            this.set_relevant_coll_balls_wall(player_dist_center);
        }
    }
    circular_move(CCW) {
        // get radians
        var pos_in_rad = coord_to_rad(this.pos);
        // get distance
        var dist = distance(this.pos, center_coord);
        // add radians
        var new_rad = pos_in_rad;
        if (CCW) { new_rad -= 0.05; } // TODO standardize on ring diameter
        if (!CCW) { new_rad += 0.05; }
        // convert back to coord
        this.pos = get_exact_coord(new_rad, dist, false);
    }
    resolve_collisions(old_pos) {

        // TODO possible problem: zero velocity but some small remaining overlap
        // due to rounding errors

        // first reset position to old pos
        // TODO: delete since not really necessary
        this.pos = {x: old_pos.x, y: old_pos.y};
        // split the velocity according to number of computations
        var vel_step = {x: this.vel.x/collision_steps, y: this.vel.y/collision_steps};
        // run the collision check collision_steps times --> only as long as no obstacle is reached
        var obstacle = false;
        // cycle through every velocity component
        for (let i = 0; i < collision_steps; i++) {
            // only update pos if no obstacle
            if (!obstacle) {
                // step forward the position
                this.pos.x += vel_step.x;
                this.pos.y += vel_step.y;

                // circles thru the collision balls
                for (let j = 0; j < this.collision_balls.length; j++) {
                    const cball = this.collision_balls[j];
                    // check if overlap --> collision needs to be resolved
                    if (distance(this.pos, cball) < (this.radius + ring_thickness/2)) {
                        // assuming player movement caused collision --> vel > 0
                        var overlap = (this.radius + ring_thickness/2) - distance(this.pos, cball);
                        // scale velocity vector with overlap and substract
                        var unit_vec_vel = unit_vector(this.vel);
                        var scaled_vel = {x: unit_vec_vel.x*overlap, y: unit_vec_vel.y*overlap};
                        this.pos.x -= 1.5*scaled_vel.x;
                        this.pos.y -= 1.5*scaled_vel.y;
                        // TODO make sure whether there are no remaining overlaps
                        // if (distance(this.pos, cball) < (this.radius + ring_thickness/2)) {
                        //     console.log("overlaps remain");
                        // }
                        // stop the collision computation steps
                        obstacle = true;
                    }
                }
            }
        }
    }
    update() {
        // copy old pos (before changes) for velocity derivation
        var old_pos = {x: this.pos.x, y: this.pos.y};

        // apply movements
        if (this.leftmove) { this.pos.x -= player_speed; }
        if (this.rightmove) { this.pos.x += player_speed; }
        if (this.upmove) { this.pos.y -= player_speed; }
        if (this.downmove) { this.pos.y += player_speed; }

        // debug: circular movements (keys Q, E)
        if (this.CWmove) { this.circular_move(false); }
        if (this.CCWmove) { this.circular_move(true); }

        // derive vel
        this.vel.x = this.pos.x - old_pos.x;
        this.vel.y = this.pos.y - old_pos.y;

        // collision resolution (in multiple steps)
        this.resolve_collisions(old_pos);

        // update current position in terms of ring indices
        this.set_current_ring_ind();

        // update virtual collision objects
        // this.set_closest_collision_balls();
        this.set_relevant_coll_balls();

        // increase in size
        // this.radius += size_grow;
    }
    render() {

        // debug: render collision_balls
        for (let index = 0; index < this.collision_balls.length; index++) {
            const cball = this.collision_balls[index];
            draw_circ(ring_thickness/2, cball, "white");
        }

        draw_circ(this.radius, this.pos, this.color);

        // debug: draw vel
        draw_line([this.pos, {x: this.pos.x + this.vel.x, y: this.pos.y + this.vel.y}], "white", 3);

    }
}

class Maze {
    // TODO replaces old rings array
    constructor(network) {
        this.network = network;
        this.rings = [];
        this.randomly_init_rings();
        this.set_coords_of_innermost_ring();
        this.resolve_ring_problems();
    }
    randomly_init_rings() {
        for (let index = 0; index < this.network.layers.length; index++) {
            var new_ring = new Ring(index, Math.random(), this.network.layers[index].length);
            this.rings.push(new_ring);
        }
    }
    set_coords_of_innermost_ring() {
        // set maximum rotation if adjacent rings have same number of holes
        // special case innermost ring: rotation does not matter
        var coords = circumference_coords_simpler(this.network.layers[0].length, start_rot);
        this.rings[0].set_holes(coords);
    }
    resolve_overlapping_holes() {

        const layer = this.network.layers[hole_problem_level];
        start_rot = (start_rot + 0.005)%1; // TODO set rotation angle proportional to circumference
        var coords = circumference_coords_simpler(layer.length, start_rot);
        this.rings[hole_problem_level].set_holes(coords);
        var this_hcoords = this.rings[hole_problem_level - 1].holes_coords;
        var inner_hcoords = this.rings[hole_problem_level].holes_coords;
        if (!this.overlap_holes(this_hcoords, inner_hcoords, hole_problem_level - 1)) {
            var relevant_id = 0;
            for (let index = 0; index < this.network.layers[hole_problem_level - 1].length; index++) {
                var parent_node = this.network.layers[hole_problem_level - 1][index];
                if (parent_node.virtual) {
                    relevant_id = parent_node.id;
                }
                if (hole_problem_level < n_rings - 1 && this.overlap_children(parent_node)) { // not last ring
                    // success
                    start_rot = Math.random();
                    // already add the walls of this layer
                    this.set_walls_one_layer(hole_problem_level);
                    this.remove_virtuals_debug_layer(hole_problem_level);
                    hole_problem_level++;
                    this.rings[parent_node.ind].overlaps = [];
                    this.rings[parent_node.ind + 1].debug_holes_inner = [];
                    break;
                }
            }
        }
    }
    set_walls_one_layer(index) {
        var tol = (this.rings[index].diameter/2)*0.00005;
        for (let index2 = 0; index2 < this.network.layers[index - 1].length; index2++) {
            const node = this.network.layers[index - 1][index2];
            // get position (simpler algorithm)
            // get correct indices
            var min_id = node.children[0].id;
            var max_id = node.children[node.children.length-1].id;
            // then convert holes coords to start/end of actual holes
            var hcoords_parent = holes_from_coords(this.rings[index - 1].holes_coords);
            var h_coords_children_layer = holes_from_coords(this.rings[index].holes_coords);
            var min = h_coords_children_layer[min_id].start;
            var max = h_coords_children_layer[max_id].end;
            // assuming parent start is always between bounds of children
            this.rings[index].add_wall((min - tol)%(2*Math.PI));
        }
    }
    remove_virtuals_debug_layer(index) {

        var non_virtual_holes = [];
        var virtual_holes = [];
    
        // convert to holes coords
        var hcoords = holes_from_coords(this.rings[index].holes_coords);
    
        // omit virtual nodes
        for (let index2 = 0; index2 < hcoords.length; index2++) {
            const node = this.network.layers[index][index2];
            if (!node.virtual) {
                non_virtual_holes.push(hcoords[index2]);
            } else {
                virtual_holes.push(hcoords[index2]);
            }
        }
    
        // debug
        this.rings[index].debug_holes = non_virtual_holes;
        this.rings[index].debug_holes_inner = virtual_holes;
    
    }
    remove_virtuals() {
        for (let index = 0; index < this.rings.length; index++) {
            // BIGTODO
            
            var non_virtual_holes = [];
            var virtual_holes = [];
    
            // convert to holes coords
            var hcoords = holes_from_coords(this.rings[index].holes_coords);
    
            // omit virtual nodes
            for (let index2 = 0; index2 < hcoords.length; index2++) {
                const node = this.network.layers[index][index2];
                if (!node.virtual) {
                    non_virtual_holes.push(hcoords[index2]);
                } else {
                    virtual_holes.push(hcoords[index2]);
                }
            }
    
            this.rings[index].set_holes_hcoords(non_virtual_holes);
        }
    }
    resolve_overlaps_last_layer(parent_node) {
        const layer = this.network.layers[hole_problem_level]; // children layer
        start_rot = (start_rot + 0.005)%1; // TODO set rotation angle proportional to circumference
        var coords = circumference_coords_simpler(layer.length, start_rot);
        this.rings[hole_problem_level].set_holes(coords);
        var this_hcoords = this.rings[hole_problem_level - 1].holes_coords;
        var inner_hcoords = this.rings[hole_problem_level].holes_coords;
        if (!this.overlap_holes(this_hcoords, inner_hcoords, hole_problem_level - 1)) {
    
            // var parent_node = layers[hole_problem_level - 1][index];
    
            if (this.overlap_children(parent_node)) {
                // success
                // already add the walls of this layer
                this.set_walls_one_layer(hole_problem_level);
                this.remove_virtuals_debug_layer(hole_problem_level);
                hole_problem_level++;
                this.rings[parent_node.ind].overlaps = [];
                this.rings[parent_node.ind + 1].debug_holes_inner = [];
            } else {
                this.rings[parent_node.ind + 1].debug_holes_inner = [this_hcoords[parent_node.id]];
            }
    
        }
    }
    overlap_holes(this_coords, new_coords, inner_layer_ind) {

        var output = false; // flags if overlaps detected
    
        var inner_holes = holes_from_coords(this_coords);
        var this_holes = holes_from_coords(new_coords);
    
        for (let index = 0; index < inner_holes.length; index++) {
            var inner_coord = inner_holes[index];
            for (let index2 = 0; index2 < this_holes.length; index2++) {
            // for (let index2 = 1; index2 < 3; index2++) {
                var this_coord = this_holes[index2];
                var problem = overlap_two_holes(this_coord, inner_coord);
    
                if (problem[0]) {
                    output = true;
                }
            }
        }
    
        return output;
    }
    overlap_children(node) {
        // 1. get rad of parent_node.start
        var parent_layer_ind = node.ind;
        var parent_id = node.id;
        var parent_hcoords = holes_from_coords(this.rings[parent_layer_ind].holes_coords);
        var rad = parent_hcoords[parent_id].start;
        // 2. get children bounds
        var h_coords_children_layer = holes_from_coords(this.rings[parent_layer_ind + 1].holes_coords);
        // 2.1 case: only one child --> set larger bounds
        if (node.children.length == 1) {
            var child_id = node.children[0].id;
            var n_nodes_child_layer = this.network.layers[parent_layer_ind + 1].length;
            var next_child_id = (child_id + 1)%n_nodes_child_layer; // handle wrap-around case
            var min = h_coords_children_layer[child_id].start;
            var max = h_coords_children_layer[next_child_id].start;
        }
        // 2.2 case: more than one child
        else {
            // get correct indices
            var min_id = node.children[0].id;
            var max_id = node.children[node.children.length-1].id;
            // then convert holes coords to start/end of actual holes
            var min = h_coords_children_layer[min_id].start;
            var max = h_coords_children_layer[max_id].end;
        }
        var bounds = {start: min, end: max};
        // 3. test if overlap
        var test = rad_between_bounds(rad, bounds);
        // debug
        if (!test) {
            this.rings[parent_layer_ind].overlaps = [rad];
            this.rings[parent_layer_ind + 1].debug_holes_inner = [bounds];
        }
        return test;
        // }
    }
    resolve_ring_problems() {
        while (hole_problem_level < n_rings) {
            if (hole_problem_level < n_rings - 1) {
                this.resolve_overlapping_holes();
            } else if (hole_problem_level == n_rings - 1) {
                // get relevant node
                var rel_ind = 0;
                for (let index = 0; index < this.network.layers[hole_problem_level - 1].length; index++) {
                    const node = this.network.layers[hole_problem_level - 1][index];
                    if (!node.virtual && !node.dead_end) {
                        rel_ind = index;
                        break;
                    }
                }
                this.resolve_overlaps_last_layer(this.network.layers[hole_problem_level - 1][rel_ind]);
            }
        }
        this.remove_virtuals();
    }
    render() {
        // maze consists of walls, rings
        // currently, walls belong to rings
        for (var i=0; i<this.rings.length; i++) {
            this.rings[i].render();
        }
    }
}

class Cell {
    constructor(ring_ind, ind_in_ring) {
        this.ring_ind = ring_ind
        this.ind_in_ring = ind_in_ring;
        // neighbours (also cells)
        this.neighbours = [];
        this.visited = false;
        this.color = "orange";
    }
    get_neighbours() {
        var outer_n = true;
        var inner_n = true;
        var this_ring = network.rings[this.ring_ind]
        var this_l = this_ring.length;
        var neighbours = [];
        // no outer neighbours if last ring
        if (this.ring_ind == n_rings - 1) { outer_n = false; }
        // no inner neighbours if first ring
        if (this.ring_ind == 0) { inner_n = false; }
        // determine outer neighbours
        if (outer_n) {
            var outer_l = network.rings[this.ring_ind + 1].length;
            // same ind if holes not doubled
            if (this_l == outer_l) {
                var outer_n_c = new Cell(this.ring_ind + 1, this.ind_in_ring);
                neighbours.push(outer_n_c);
            }
            // else: doubled - two outer neighbours
            else {
                var outer_n_c1 = new Cell(this.ring_ind + 1, this.ind_in_ring*2);
                var outer_n_c2 = new Cell(this.ring_ind + 1, this.ind_in_ring*2 + 1);
                neighbours.push(outer_n_c1);
                neighbours.push(outer_n_c2);
            }
        }
        // determine inner neighbours
        if (inner_n) {
            var inner_l = network.rings[this.ring_ind - 1].length;
            // - same ind if holes not doubled
            if (this_l == inner_l) {
                var inner_n_c = new Cell(this.ring_ind - 1, this.ind_in_ring);
                neighbours.push(inner_n_c);
            }
            // else: halved - different ind
            else {
                var inner_n_c = new Cell(this.ring_ind - 1, Math.floor(this.ind_in_ring/2));
                neighbours.push(inner_n_c);
            }
        }
        // determine sidewards neighbours
        if (this_l == 2) {
            // only left
            if (this.ind_in_ring + 1 == this_l) {
                var s_n = new Cell(this.ring_ind, 0);
            } else {
                var s_n = new Cell(this.ring_ind, this.ind_in_ring + 1);
            }
            neighbours.push(s_n);
        } else if (this_l > 2) {
            // left
            if (this.ind_in_ring + 1 == this_l) {
                var s_n1 = new Cell(this.ring_ind, 0);
            } else {
                var s_n1 = new Cell(this.ring_ind, this.ind_in_ring + 1);
            }
            // right
            if (this.ind_in_ring - 1 < 0) {
                var s_n2 = new Cell(this.ring_ind, this_l - 1);
            } else {
                var s_n2 = new Cell(this.ring_ind, this.ind_in_ring - 1);
            }
            neighbours.push(s_n1);
            neighbours.push(s_n2);
        }
        return neighbours;
    }
    get_rad() {
        var rads = equidistant_rads(network.rings[this.ring_ind].length);
        return rads[this.ind_in_ring];
    }
    get_pos() {
        // transform into coordinate
        // get distance from ring_ind
        var dist = this.ring_ind*ring_spacing + 0.5*ring_spacing;
        var rad = this.get_rad();
        return get_exact_coord(rad, dist, false);
    }
    render() {
        draw_circ(5, this.get_pos(), this.color);
    }
}

class MazeNetworkRandom {
    constructor() {
        // certain number of cells in each ring
        this.rings = [];
        // keeping track if number of holes should get doubled
        this.current_holes_n = 0;
        this.n_cells = 0;
        this.init_cells();
        this.current_cell = this.rings[0][0];
        this.previous_cell = this.rings[0][0];
        this.visited_cells = []; // stores the occupied cells
        this.stack = [];
        this.neighbour_cells = [];
        this.lines = [];
    }
    init_cells() {
        // BIG TODO: adjust so that n_cells == n_walls, always!
        // exception first ring: 4 cells even though no walls
        var n_holes = this.get_number_of_holes_ring(0);
        var cells = [];
        for (let index2 = 0; index2 < n_holes; index2++) {
            var cell = new Cell(0, index2);
            cells.push(cell);
            this.n_cells++;
        }
        this.rings.push(cells);
        for (let index = 1; index < n_rings; index++) {
            var cells = [];
            // get number of holes and set whether doubled
            var n_holes = this.get_number_of_holes_ring(index - 1);
            for (let index2 = 0; index2 < n_holes; index2++) {
                var cell = new Cell(index, index2);
                cells.push(cell);
                this.n_cells++;
            }
            this.rings.push(cells);
        }
    }
    set_closest_cell(pos) {
        // for debugging
        var current_dist = Infinity;
        var current_closest_cell = this.rings[0][0];
        for (let index = 0; index < this.rings.length; index++) {
            const ring = this.rings[index];
            for (let index2 = 0; index2 < ring.length; index2++) {
                const cell_pos = ring[index2].get_pos();
                if (distance(cell_pos, pos) < current_dist) {
                    current_dist = distance(cell_pos, pos);
                    current_closest_cell = ring[index2];
                }
            }
        }
        this.previous_cell = this.current_cell; // TODO: make copy?
        this.current_cell = current_closest_cell;
    }
    get_number_of_holes_ring(ring_ind) {
        // TODO some problem here - too many holes detected sometimes
        // how many holes (length r) fit onto circumference?
        var r = new Ring(ring_ind);
        var c = circumference_from_diameter(r.diameter);
        var n = Math.floor(c/ring_spacing);
        // make sure the walls fit as well
        while (n*ring_spacing + n*ring_thickness >= c) {
            n--;
        }
        // set to current holes
        if (this.current_holes_n == 0) {
            this.current_holes_n = n;
        } else if (n >= this.current_holes_n*2) {
            this.current_holes_n = this.current_holes_n*2;
        }
        return this.current_holes_n;
    }
    render() {
        // debug: render visited cells
        // check if last line already appended
        if (this.lines.length > 0) {
            var lastline = this.lines[this.lines.length-1];
            if (lastline[1].x != this.current_cell.get_pos().x || 
                lastline[1].y != this.current_cell.get_pos().y) {
                    var pos1 = this.previous_cell.get_pos();
                    var pos2 = this.current_cell.get_pos();
                    this.lines.push([pos1, pos2]);
            }
        } else {
            var pos1 = this.previous_cell.get_pos();
            var pos2 = this.current_cell.get_pos();
            this.lines.push([pos1, pos2]);
        }
        for (let index = 0; index < this.visited_cells.length; index++) {
            const visited_cell = this.visited_cells[index];
            visited_cell.color = "red";
            visited_cell.render();
            // draw line to current cell as well
            // if (index == this.visited_cells.length - 1) {
            //     var pos1 = visited_cell.get_pos();
            //     var pos2 = this.current_cell.get_pos();
            //     draw_line([pos1, pos2], "black");
            // }
        }
        for (let index = 0; index < this.lines.length; index++) {
            draw_line(this.lines[index], "black");
        }
        // debug: render current neighbour cells (set in maze.dig function)
        for (let index = 0; index < this.neighbour_cells.length; index++) {
            const neighbour = this.neighbour_cells[index];
            neighbour.color = "white";
            neighbour.render();
        }
        this.current_cell.color = "orange";
        this.current_cell.render();
    }
}

class MazeRandom {
    constructor(network) {
        // get info from network
        this.network = network;
        // generate the rings
        this.rings = [];
        this.go_back = false; // debug
        this.set_rings();
        this.init_walls();
        this.dig_all();
    }
    set_rings() {
        // already sets the holes accordingly
        for (let index = 0; index < n_rings; index++) {
            // add ring to array
            var r = new Ring(index);
            var c = circumference_from_diameter(r.diameter);
            // get number of holes
            var n = this.network.rings[index].length;
            // set holes
            // var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
            // hcoords = [];
            // r.set_holes_hcoords(hcoords);
            // add to rings array
            this.rings.push(r);
        }
    }
    init_walls() {
        // TODO: subdivision should occur at higher level
        for (let index = 1; index < this.rings.length; index++) {
            var n_walls = this.network.rings[index].length;
            var wall_coords = equidistant_wall_coords(n_walls);
            for (let index2 = 0; index2 < wall_coords.length; index2++) {
                var wall = new Wall(index - 1, index, wall_coords[index2]);
                this.rings[index].walls.push(wall);
            }
        }
    }
    remove_wall(ring_ind, wall_ind) {
        this.rings[ring_ind].walls[wall_ind].dug = true;
    }
    exclude_visited_neighbours(neighbours) {
        // TODO figure out why this returns one too many in case of going back
        var reduced_neighbours = [];
        for (let index = 0; index < neighbours.length; index++) {
            const nb = neighbours[index];
            var valid = true;
            for (let index2 = 0; index2 < network.visited_cells.length; index2++) {
                const vc = network.visited_cells[index2];
                if (nb.ring_ind == vc.ring_ind && nb.ind_in_ring == vc.ind_in_ring) {
                    valid = false;
                    break;
                }
            }
            if (valid) {
                reduced_neighbours.push(nb);
            }
        }
        return reduced_neighbours;
    }
    set_neighbours() {
        // first randomly select neighbouring cell
        var nbs = this.network.current_cell.get_neighbours();
        // make sure visited cells are excluded
        var nbs_MINUS_v = this.exclude_visited_neighbours(nbs);
        // debug: set the neighbour cells as property to network
        this.network.neighbour_cells = nbs_MINUS_v;
        return nbs_MINUS_v;
    }
    set_current_cell_to_last_valid_pos() {
        // goes through stack
        // determines number of valid neighbours
        // sets current cell as soon as there are >0 valid neighbours
        // also removes from stack
        // pop last element from stack
        network.stack.pop();
        const visited_cell = network.stack[network.stack.length - 1];
        network.previous_cell = network.current_cell; // TODO: make copy?
        network.current_cell = visited_cell;
        // determine neighbours
        var n_n = this.set_neighbours().length;
        if (n_n > 0) {
            // TODO figure out why 1 must be removed
            // stop search
            this.go_back = false;
            this.dig();
        }
    }
    set_visited_cell(c) {
        // first special case: no cells visited (start of program)
        if (network.stack.length == 0) {
            // set the current cell as visited
            c.visited = true;
            network.visited_cells.push(c);
            network.stack.push(c);
        } else {
            var last_visited = network.stack[network.stack.length - 1];
            if (c.ring_ind != last_visited.ring_ind || c.ind_in_ring != last_visited.ind_in_ring) {
                // set the current cell as visited no matter what
                c.visited = true;
                network.visited_cells.push(c);
                network.stack.push(c);
            }
        }
    }
    dig() {
        // removes walls or creates holes
        // TODO: put in correct place (debug)
        if (this.go_back) {
            this.set_current_cell_to_last_valid_pos();
            return;
        }
        // set cell as visited if not already set
        var c = network.current_cell;
        this.set_visited_cell(c);
        var n = this.set_neighbours();
        if (n.length == 0) {
            // set current cell to last location with options
            // but still store this cell as visited in list
            this.go_back = true;
            return;
        }
        // select randomly from cells
        var r_i = Math.floor(Math.random() * Math.floor(n.length));
        var rand_i = n[r_i];
        // introduce abbreviations
        var rri = rand_i.ring_ind;
        var cri = network.current_cell.ring_ind;
        var rir = rand_i.ind_in_ring;
        var cir = network.current_cell.ind_in_ring;
        var n_cells = network.rings[cri].length;
        var n_walls = this.rings[cri].walls.length;
        // sidewards: remove wall
        if (rri == cri && rri > 0) {
            if (rir == 0 && cir == network.rings[cri].length - 1) {
                // CASE 1
                // new cell right of old cell but wrapped around
                var dir = 1;
                var rwi = 0;
            } else if (rir == network.rings[cri].length - 1 && cir == 0) {
                // CASE 2
                // new cell left of old cell but wrapped around
                var dir = -1;
                var rwi = 0;
            } else if (rir < cir) {
                // CASE 3
                // new cell left of old cell and not wrapped around
                var dir = -1;
                var rwi = rir + 1;
            } else if (rir > cir) {
                // CASE 4
                // new cell right of old cell and not wrapped around
                var dir = 1;
                var rwi = rir;
            }
            // handle case: n_cells == 2*n_walls
            if (n_cells > n_walls) { // doubled
                if (dir == -1) {
                    if (rir == n_cells - 1) {
                        rwi = 0;
                    } else if (rir%2 == 1) {
                        rwi = Math.ceil(rir/2);
                    } else {
                        rwi = -1;
                    }
                } else {
                    if (rir%2 == 0) {
                        rwi = Math.floor(rir/2);
                    } else {
                        rwi = -1;
                    }
                }
            }
            // only set as dug if wall is to be removed
            if (rwi >= 0) {
                this.rings[rri].walls[rwi].dug = true;
            }
        }
        // through layers: add holes
        if (rri != cri) {
            // identify which ring
            var ri = Math.min(rri, cri);
            // setup holes coords to choose from
            var ring = network.rings[ri + 1];
            var n = ring.length;
            var c = this.rings[ri].get_c();
            var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
            // identify ind in ring
            if (rri < cri) {
                var ir = cir;
            } else {
                var ir = rir;
            }
            this.rings[ri].add_hcoord_to_hcoords(hcoords[ir]);
        }
        // set current cell to new position
        network.previous_cell = network.current_cell; // TODO: make copy?
        network.current_cell = rand_i;
    }
    debug_dig(dir) {
        // debug: dir: -1 = left, 1 = right
        // stay in same ring space --> only remove walls
        // TODO only allow this function to 
        var c = network.current_cell;
        var cri = network.current_cell.ring_ind;
        var cir = network.current_cell.ind_in_ring;
        // TODO: consider case when n_walls < n_cells
        var n_cells = network.rings[cri].length;
        var n_walls = this.rings[cri].walls.length;
        if (dir < 0 && cir == 0) {
             // CASE 2
            var wall_rmv_ind = 0;
            var rir = network.rings[cri].length - 1;
        } else if (dir > 0 && cir == network.rings[cri].length - 1) {
            // CASE 1
            var wall_rmv_ind = 0;
            var rir = 0;
        } else if (dir < 0) {
            var wall_rmv_ind = (cir + dir + 1)%network.rings[cri].length;
            var rir = cir + dir;
        } else {
            var wall_rmv_ind = cir + dir;
            var rir = cir + dir;
        }
        if (n_cells > n_walls) { // doubled
            if (dir == -1) {
                if (rir == n_cells - 1) {
                    wall_rmv_ind = 0;
                } else if (rir%2 == 1) {
                    wall_rmv_ind = Math.ceil(rir/2);
                } else {
                    wall_rmv_ind = -1;
                }
            } else {
                if (rir%2 == 0) {
                    wall_rmv_ind = Math.floor(rir/2);
                } else {
                    wall_rmv_ind = -1;
                }
            }
        }
        if (wall_rmv_ind >= 0 && cri > 0) {
            this.rings[cri].walls[wall_rmv_ind].dug = true;
        }
        // store previous cell and set current cell to new position
        c.visited = true;
        network.visited_cells.push(c);
        network.previous_cell = network.current_cell; // TODO: make copy?
        network.current_cell = new Cell(cri, rir);
    }
    dig_all() {
        var all_dug = false;
        while (!all_dug) {
            this.dig();
            if (network.visited_cells.length == network.n_cells) {
                all_dug = true;
            }
        }
        // set the last (escape) hole
        var ring = network.rings[network.rings.length-1];
        var n = ring.length;
        var c = this.rings[this.rings.length-1].get_c();
        var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
        var ir = Math.floor(Math.random()*hcoords.length);
        console.log(ir);
        this.rings[this.rings.length-1].add_hcoord_to_hcoords(hcoords[ir]);
    }
    update() {

        // if (network.visited_cells.length == network.n_cells) {
        //     console.log("ALL DUG");
        //     // set the last (escape) hole
        //     var ring = network.rings[network.rings.length-1];
        //     var n = ring.length;
        //     var c = this.rings[this.rings.length-1].get_c();
        //     var hcoords = equidistant_holes_coords(n, c, ring_spacing, ring_thickness);
        //     var ir = Math.random()*hcoords.length;
        //     this.rings[this.rings.length-1].add_hcoord_to_hcoords(hcoords[ir]);
        // } else {
        //     this.dig();
        // }
    }
    render() {
        // maze consists of walls, rings
        // currently, walls belong to rings
        for (var i=0; i<this.rings.length; i++) {
            this.rings[i].render();
        }
    }
}

class Wall {
    constructor(ind1, ind2, rad, invalid=false, ring1, ring2) {
        this.ind1 = ind1; // inner ring to which wall is connected
        this.ind2 = ind2; // this ring to which wall is connected
        this.ring1 = ring1;
        this.ring2 = ring2;
        this.rad = rad; // angle in radians
        // derived vars
        this.r1 = ((ind1+1)*2*ring_spacing - ring_thickness)/2;
        this.r2 = ((ind2+1)*2*ring_spacing - ring_thickness)/2;
        this.coord1 = get_exact_coord(this.rad, this.r1, false);
        this.coord2 = get_exact_coord(this.rad, this.r2, false);
        this.invalid = invalid; // for debugging
        this.moving = false;
        this.dug = false;
    }
    set_validity() {
        // get them coords
        var coords1 = this.ring1.holes_coords;
        var coords2 = this.ring2.holes_coords;
        // if any overlap with hole of either inner or outer ring --> function sets validity
        this.invalid = false;
        if (overlap_hole_region(this.rad, coords1)) { this.invalid = true; }
        if (overlap_hole_region(this.rad, coords2)) { this.invalid = true; }
    }
    update_coords() {
        this.coord1 = get_exact_coord(this.rad, this.r1, false);
        this.coord2 = get_exact_coord(this.rad, this.r2, false);
    }
    update() {
        this.set_validity();
        // position updates (debugging)
        if (this.moving) {
            this.rad -= 0.001;
            this.rad %= Math.PI*2; // make sure value stays within bounds
            this.update_coords();
        }
    }
    render() {
        var color = "black"
        if (!this.dug) {
            draw_line([this.coord1, this.coord2], color, ring_thickness);
        }
    }
}

class Ring {
    constructor(ind) {
        // lower inds --> smaller rings
        this.ind = ind;
        this.diameter = (ind+1)*2*ring_spacing - ring_thickness;
        this.holes_coords = [{start: 0, end: Math.PI*2}]; // init without holes
        this.edge_ball_coords = [];
        this.walls = [];
        this.rotated = false; // flag to store whether holes are rotated
        this.overlaps = [];
        this.overlaps_inner = [];
        this.debug_holes = [];
        this.debug_holes_inner = [];
        this.holes_overlap = false;
    }
    get_c() {
        // circumference from diameter
        return Math.PI*this.diameter;
    }
    set_holes(coords) {
        var coords_container = circumference_coords_given_coords(coords, this.diameter/2, ring_spacing, ring_thickness);
        this.holes_coords = coords_container[0];
        this.edge_ball_coords = coords_container[1];
    }
    add_coord_to_hcoords(coord) {
        // first remove placeholder:
        if (this.holes_coords.length == 1) {
            if (this.holes_coords[0].start == 0 && this.holes_coords[0].end == Math.PI*2) {
                this.holes_coords = [];
            }
        }
        var c = this.diameter*Math.PI;
        // flip coord (not sure why necessary)
        // coord = flip_rad(coord);
        var new_hcoord = start_end_from_coord(coord, c, ring_spacing, ring_thickness, true);
        this.holes_coords.push(new_hcoord);
        this.edge_ball_coords.push(new_hcoord.start);
        this.edge_ball_coords.push(new_hcoord.end);
    }
    add_hcoord_to_hcoords(hcoord) {
        // first remove placeholder:
        if (this.holes_coords.length == 1) {
            if (this.holes_coords[0].start == 0 && this.holes_coords[0].end == Math.PI*2) {
                this.holes_coords = [];
            }
        }
        this.holes_coords.push(hcoord);
        console.log(hcoord.start);
        // convert to coords
        var start = get_exact_coord(hcoord.start, this.diameter/2, false);
        var end = get_exact_coord(hcoord.end, this.diameter/2, false);
        this.edge_ball_coords.push(start);
        this.edge_ball_coords.push(end);
    }
    set_holes_hcoords(hcoords) {
        this.holes_coords = [];
        this.edge_ball_coords = [];
        for (let index = 0; index < hcoords.length; index++) {
            this.edge_ball_coords.push(hcoords[index].start);
            this.edge_ball_coords.push(hcoords[index].end);
        }
        this.holes_coords.push({start: hcoords[hcoords.length-1].end, end: hcoords[0].start});
        for (let index = 0; index < hcoords.length - 1; index++) {
            this.holes_coords.push({start: hcoords[index].end, end: hcoords[index + 1].start});
        }
    }
    /*
    set_holes_hcoords2() {
        if (this.holes_coords.length == 0) {
            this.holes_coords.push({start: 0, end: Math.PI*2});
        } else {
            this.edge_ball_coords = [];
            for (let index = 0; index < this.holes_coords.length; index++) {
                this.edge_ball_coords.push(this.holes_coords[index].start);
                this.edge_ball_coords.push(this.holes_coords[index].end);
            }
            this.holes_coords.push({start: this.holes_coords[hcoords.length-1].end, end: this.holes_coords[0].start});
            for (let index = 0; index < hcoords.length - 1; index++) {
                this.holes_coords.push({start: hcoords[index].end, end: hcoords[index + 1].start});
            }
        }
    }
    */
    add_wall(rad) {
        var lower_level_ring = rings[this.ind - 1];
        this.walls.push(new Wall(this.ind, this.ind-1, rad, true, this, lower_level_ring));
    }
    add_walls(lower_level_ring) {
        // /*
        // repeat finding random wall radians position until no overlap with holes
        var coords1 = this.holes_coords;
        var coords2 = lower_level_ring.holes_coords;
        var wall_rads = random_wall_rads(Math.ceil((this.ind + 1)*1.3));
        // var wall_rads = random_wall_rads(1);
        var validities = get_wall_validity(coords1, coords2, wall_rads);
        for (let index = 0; index < wall_rads.length; index++) {
            this.walls.push(new Wall(this.ind, this.ind-1, wall_rads[index], validities[index], this, lower_level_ring));
        }
        // */
    }
    update() {
        for (let index = 0; index < this.walls.length; index++) {
            this.walls[index].update();
        }
        for (let index = 0; index < this.walls.length; index++) {
            if (this.walls[index].invalid) {
                this.walls[index].moving = true;
                this.walls[index].update();
            } else {
                this.walls[index].moving = false;
            }
        }
    }
    render() {
        // convert holes_coords to line segment coords
        var line_seg_coords = coords_from_holes(this.holes_coords);
        // draw line segments
        for (let index = 0; index < line_seg_coords.length; index++) {
            draw_circ_segment(this.diameter/2, {x: canvas.width/2, y: canvas.height/2}, "black", ring_thickness, line_seg_coords[index]);
        }
        // draw edges
        for (let index = 0; index < this.edge_ball_coords.length; index++) {
            draw_circ(ring_thickness/2, this.edge_ball_coords[index], "black");
        }
        // draw walls
        for (let index = 0; index < this.walls.length; index++) {
            this.walls[index].render();
        }
    }
}

// instantiate objects
add_values();
var player = new Player();
var network = new MazeNetworkRandom();
var maze = new MazeRandom(network);

// overall update function
function update() {
    // run changes (update objects)
    player.update();
    // for (let index = 0; index < rings.length; index++) {
    //     rings[index].update();
    // }

    // draw all changes
    draw();
    // get animation going
    requestAnimationFrame(update);
}

// overall draw function
function draw() {
    // refresh
    set_canvas_bg("lightblue");
    // draw network (currently debug)
    if (debug_draw_toggle) {
        network.render();
    }
    // draw maze
    maze.render();
    // draw player
    player.render();
}

// event listener actions
function keydown(e) {
    // left
    if (e.keyCode == 37) { player.leftmove = true; }
    // up
    if (e.keyCode == 38) { player.upmove = true; }
    // right
    if (e.keyCode == 39) { player.rightmove = true; }
    // down
    if (e.keyCode == 40) { player.downmove = true; }

    // enter --> debug rotate walls
    if (e.keyCode == 13) {
        maze.update();
    }

    // S, W --> debug increase/decrease player size
    if (e.keyCode == 87) { player.radius += 20*size_grow; }
    if (e.keyCode == 83 && player.radius > 1 + size_grow) { player.radius -= 20*size_grow; }
    // Q, E --> make circular movements
    if (e.keyCode == 81) { player.CCWmove = true; }
    if (e.keyCode == 69) { player.CWmove = true; }
    // A, D --> debug erase walls
    if (e.keyCode == 65) { maze.debug_dig(-1); }
    if (e.keyCode == 68) { maze.debug_dig(1); }
    // B --> toggle debug
    if (e.keyCode == 66) { debug_draw_toggle = !debug_draw_toggle; }

}
function keyup(e) {
    // left
    if (e.keyCode == 37) { player.leftmove = false; }
    // up
    if (e.keyCode == 38) { player.upmove = false; }
    // right
    if (e.keyCode == 39) { player.rightmove = false; }
    // down
    if (e.keyCode == 40) { player.downmove = false; }

    // enter
    if (e.keyCode == 13) { 
        for (let index = 0; index < rings.length; index++) {
            for (let index2 = 0; index2 < rings[index].walls.length; index2++) {
                rings[index].walls[index2].moving = false;
            }
        }
    }

    // Q, E --> make circular movements
    if (e.keyCode == 81) { player.CCWmove = false; }
    if (e.keyCode == 69) { player.CWmove = false; }
}
function mouseclick(e) {
    // debug: select cell
    var pos = getXY(e);
    player.pos.x = pos.x;
    player.pos.y = pos.y;
}

// start game loop
update();