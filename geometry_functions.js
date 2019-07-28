var center_coord = {x: canvas.width/2, y: canvas.height/2};

function getXY(e) {
	
	// return exact position relative to canvas (for physics)
	var x_exact = mouse_pos_relative_to_canvas(e).x;
	var y_exact = mouse_pos_relative_to_canvas(e).y;
	return {x: x_exact, y: y_exact};
	
}

function mouse_pos_relative_to_canvas(e) {
	
	var xval = e.clientX - canvas.offsetLeft;
	var yval = e.clientY - canvas.offsetTop;
	return {x: xval, y: yval};
	
}

function deg_to_rad(deg) {
    return deg*(Math.PI/180);
}

function flip_rad(rad) {
    return (rad + Math.PI) % (2*Math.PI);
}

function coords_diff(coord1, coord2) {
    var x1 = coord1.x;
    var x2 = coord2.x;
    var y1 = coord1.y;
    var y2 = coord2.y;
    var diff_x = x2 - x1;
    var diff_y = y2 - y1;
    return {x: diff_x, y: diff_y};
}

function distance(pos1,pos2) {
    var diff = coords_diff(pos1, pos2);
    // Pythagoras' theorem
    return Math.sqrt(diff.x*diff.x + diff.y*diff.y);
}

function coord_to_rad(coord) {
    // relative to center!
    var diff = coords_diff(coord, center_coord);
    return Math.atan2(diff.y, diff.x);
}

// trigonometric functions in JS are radians-based
function get_exact_coord(angle, dist, in_deg) {
    if (in_deg) {
        // TODO: flip angle
        var rad = deg_to_rad(angle); // convert
    } else {
        // flip angle (TODO: figure out why necessary)
        var rad = angle; // flip_rad(angle);
    }
    var x_val = dist*Math.cos(rad) + canvas.width/2;
    var y_val = dist*Math.sin(rad) + canvas.height/2;
    return {x: x_val, y: y_val};
}

function unit_vector(vec) {
    // get length
    var length = distance(vec, {x: 0, y: 0});
    // divide by length if not zero
    if (length > 0) {
        return {x: vec.x/length, y: vec.y/length};
    } else {
        return {x: 0, y: 0};
    }
}

function shuffle(a) {
    var j, x, i;
    for (i = a.length - 1; i > 0; i--) {
        j = Math.floor(Math.random() * (i + 1));
        x = a[i];
        a[i] = a[j];
        a[j] = x;
    }
    return a;
}

function flip_coords(coords) {
    var output = [];
    for (let index = 0; index < coords.length; index++) {
        const c = coords[index];
        output.push({start: flip_rad(c.start), end: flip_rad(c.end)});
    }
    return output;
}

function holes_from_coords(coords, flip = false) {

    var output = [];

    if (flip) {
        coords = flip_coords(coords);
    }

    // all holes but last one
    for (let j = 0; j < coords.length-1; j++) {
        var hole_start = coords[j].end %(Math.PI*2);
        var hole_end = coords[j+1].start %(Math.PI*2);
        output.push({start: hole_start, end: hole_end});
    }

    // last hole
    var hole_start = coords[coords.length-1].end %(Math.PI*2);
    var hole_end = coords[0].start %(Math.PI*2);
    output.push({start: hole_start, end: hole_end});
    
    return output;

}

function coords_from_holes(hcoords, flip=false, log=false) {
    // sort the holes in ascending manner
    var holes = sort_hcoords(hcoords, log);
    // var holes = hcoords;
    // debug: check if working
    if (log) {
        console.log(holes);
    }
    // init with wrap around coord
    var output = [{start: holes[holes.length-1].end, end: holes[0].start}];

    if (flip) {
        coords = flip_coords(coords);
    }

    // reverse process from holes_from_coords
    for (let index = 1; index < holes.length; index++) {
        output.push({start: holes[index - 1].end, end: holes[index].start});
    }

    return output;
}

function sort_hcoords(hcoords, log=false) {
    var sortable = [];
    for (let index = 0; index < hcoords.length; index++) {
    // for (var coord in hcoords) {
        sortable.push([hcoords[index].start, hcoords[index].end]);
    }
    if (log) {
        console.log(sortable);
    }
    sortable.sort(function(a, b) {
        return a[0]-b[0];
    });
    new_hcoords = [];
    for (let index = 0; index < sortable.length; index++) {
    // for (var sorted in sortable) {
        new_hcoords.push({start: sortable[index][0], end: sortable[index][1]});
    }
    return new_hcoords;
}

function get_simpler_coords(coords) {
    // i.e. only starts
    var output = [];
    for (let index = 0; index < coords.length; index++) {
        const c = coords[index];
        output.push(c.start);
    }
    return output;
}

function random_wall_rads(n_walls) {
    // returns position along ring in radians for n walls (equally spaced)
    var wall_rads = [];
    var step = (2*Math.PI)/n_walls;
    var start_rad = Math.random()*(2*Math.PI);
    for (let index = 0; index < n_walls; index++) {
        wall_rads.push(start_rad);
        start_rad += step;
        start_rad = start_rad%(2*Math.PI);
    }
    return wall_rads;
}

function overlap_hole_region(wall_rad, coords, return_coords = false, flip = true, flipc = false) {

    // flip wall rad again TODO find out why
    if (flip) {
        wall_rad = flip_rad(wall_rad);
    }
    if (flipc) {
        coords = flip_coords(coords);
    }

    // stores bool and coords
    var overlaps_hole = false;
    var hole_coords = {start: 0, end: 0};

    for (let j = 0; j < coords.length; j++) {
        var hole_start = coords[j].start;
        var hole_end = coords[j].end;
        // consider case that hole is around radians 0
        if (hole_end < hole_start) {
            // test instead if in non-hole area, then reverse return bool
            if (wall_rad >= hole_end && wall_rad <= hole_start) {
                // return false;
            } else {
                overlaps_hole = true;
                hole_coords = {start: hole_start, end: hole_end};
                break;
            }
        }
        // normal case
        if (wall_rad >= hole_start && wall_rad <= hole_end) {
            overlaps_hole = true;
            hole_coords = {start: hole_start, end: hole_end};
            break;
        }
    }

    // returns
    if (return_coords) {
        return [overlaps_hole, hole_coords];
    } else {
        return overlaps_hole;
    }
}

function overlap_hole_region_copy(wall_rad, coords, return_coords = false, flip = true, flipc = false) {

    // flip wall rad again TODO find out why
    if (flip) {
        wall_rad = flip_rad(wall_rad);
    }
    if (flipc) {
        coords = flip_coords(coords);
    }

    // stores bool and coords
    var overlaps_hole = true;
    var hole_coords = {start: 0, end: 0};

    for (let j = 0; j < coords.length-1; j++) {
        var hole_start = coords[j].end;
        var hole_end = coords[j+1].start;
        // consider case that hole is around radians 0
        if (hole_end < hole_start) {
            // test instead if in non-hole area, then reverse return bool
            if (wall_rad >= hole_end && wall_rad <= hole_start) {
                // return false;
            } else {
                overlaps_hole = false;
                hole_coords = {start: hole_start, end: hole_end};
                break;
            }
        }
        if (wall_rad >= hole_start && wall_rad <= hole_end) {
            overlaps_hole = false;
            hole_coords = {start: hole_start, end: hole_end};
            break;
        }
    }
    // test last hole (if not already a break reached)
    if (coords.length > 0 && !overlaps_hole) {
        var hole_start = coords[coords.length-1].end;
        var hole_end = coords[0].start;
        // consider case that hole is around radians 0
        if (hole_end < hole_start) {
            // test instead if in non-hole area, then reverse return bool
            if (wall_rad >= hole_end && wall_rad <= hole_start) {
                // return false;
            } else {
                overlaps_hole = false;
                hole_coords = {start: hole_start, end: hole_end};
            }
        }
        if (wall_rad >= hole_start && wall_rad <= hole_end) {
            overlaps_hole = false;
            hole_coords = {start: hole_start, end: hole_end};
        }
    }

    // returns
    if (return_coords) {
        return [overlaps_hole, hole_coords];
    } else {
        return overlaps_hole;
    }
}

function overlaps_hole(rad, coords, return_coords) {
    for (let index = 0; index < coords.length; index++) {
        
    }
}

function overlap_hole_region_OLD(wall_rad, coords, return_coords = false, flip = true, flipc = false) {

    // flip wall rad again TODO find out why
    if (flip) {
        wall_rad = flip_rad(wall_rad);
    }
    if (flipc) {
        coords = flip_coords(coords);
    }

    // stores bool and coords
    var overlaps_hole = false;
    var hole_coords = {start: 0, end: 0};

    for (let j = 0; j < coords.length-1; j++) {
        var hole_start = coords[j].end;
        var hole_end = coords[j+1].start;
        // consider case that hole is around radians 0
        if (hole_end < hole_start) {
            // test instead if in non-hole area, then reverse return bool
            if (wall_rad >= hole_end && wall_rad <= hole_start) {
                hole_coords = {start: hole_start, end: hole_end};
                // return false;
            } else {
                overlaps_hole = true;
                hole_coords = {start: hole_start, end: hole_end};
                break;
            }
        }
        if (wall_rad >= hole_start && wall_rad <= hole_end) {
            overlaps_hole = true;
            hole_coords = {start: hole_start, end: hole_end};
            break;
        }
    }
    // test last hole (if not already a break reached)
    if (coords.length > 0 && !overlaps_hole) {
        var hole_start = coords[coords.length-1].end;
        var hole_end = coords[0].start;
        // consider case that hole is around radians 0
        if (hole_end < hole_start) {
            // test instead if in non-hole area, then reverse return bool
            if (wall_rad >= hole_end && wall_rad <= hole_start) {
                hole_coords = {start: hole_start, end: hole_end};
                // return false;
            } else {
                overlaps_hole = true;
                hole_coords = {start: hole_start, end: hole_end};
            }
        }
        if (wall_rad >= hole_start && wall_rad <= hole_end) {
            overlaps_hole = true;
            hole_coords = {start: hole_start, end: hole_end};
        }
    }

    // returns
    if (return_coords) {
        return [overlaps_hole, hole_coords];
    } else {
        return overlaps_hole;
    }
}

function rad_between_bounds(rad, coords, return_coords=false) {
    // store the coord
    var overlap = false;
    // make copy
    var bounds = {start: coords.start, end: coords.end};
    // first handle case when bounds wrap around
    if (bounds.end < bounds.start) {
        // already determined when rad is even smaller --> must be within bounds
        if (rad <= bounds.end) {
            overlap = true;
        }
        // else handle case
        bounds.end += 2*Math.PI;
    }
    if (bounds.start <= rad && rad <= bounds.end) {
        overlap = true;
    }
    // any other case: not between
    if (return_coords) {
        return [overlap, bounds];
    }
    return overlap;
}

function overlap_holecoords(wall_rad, hcoords, return_coords = false, flipa = [true, false]) {

    // coordinates are not for start-end of walls but start-end of holes already

    // flip wall rad again TODO find out why
    if (flipa[0]) {
        wall_rad = flip_rad(wall_rad);
    }
    if (flipa[1]) {
        hcoords = flip_coords(hcoords);
    }

    // stores bool and coords
    var overlaps_hole = false;

    var hole_coords = {start: 0, end: 0};

    for (let j = 0; j < hcoords.length; j++) {
        // consider case that hole is around radians 0
        if (hcoords.end < hcoords.start) {
            // test instead if in non-hole area, then reverse return bool
            if (wall_rad >= hcoords.end && wall_rad <= hcoords.start) {
                // return false;
            } else {
                overlaps_hole = true;
                hole_coords = hcoords[j];
                break;
            }
        }
        if (wall_rad >= hcoords.start && wall_rad <= hcoords.end) {
            overlaps_hole = true;
            hole_coords = hcoords[j];
            break;
        }
    }

    // returns
    if (return_coords) {
        return [overlaps_hole, hole_coords];
    } else {
        return overlaps_hole;
    }
}

function overlap_two_holes(hole_coords1, hole_coords2) {
    // call the basic function four times, if any overlap --> return
    // always returns

    function find_problems(rad, coords) {
        var variant = rad_between_bounds(rad, coords, return_coords = true);
        return variant
    }

    // /*

    // INNER as thresholds
    // variant 1
    var wall_rad = hole_coords1.start;
    var coords = hole_coords2;
    var variant = find_problems(wall_rad, coords);
    if (variant[0]) { return variant; }
    // variant 2
    var wall_rad = hole_coords1.end;
    var coords = hole_coords2;
    var variant = find_problems(wall_rad, coords);
    if (variant[0]) { return variant; }

    // */

    /*

    // OUTER as thresholds
    // variant 1
    var wall_rad = hole_coords2.start;
    var coords = hole_coords1;
    var variant = find_problems(wall_rad, coords);
    if (variant[0]) { return variant; }
    // variant 2
    var wall_rad = hole_coords2.end;
    var coords = hole_coords1;
    var variant = find_problems(wall_rad, coords);
    if (variant[0]) { return variant; }

    // */

    // no problem until now
    return [false, {start: 0, end: 0}];
}

function walls_in_holes(coords1, coords2, wall_rads) {
    for (let i = 0; i < wall_rads.length; i++) {
        const wall_rad = wall_rads[i];
        // if any overlap with hole of either inner or outer ring --> function returns true
        if (overlap_hole_region(wall_rad, coords1)) { return true; }
        if (overlap_hole_region(wall_rad, coords2)) { return true; }
    }
    // no overlap
    return false;
}

function get_wall_validity(coords1, coords2, wall_rads) {
    var wall_list = [];
    for (let i = 0; i < wall_rads.length; i++) {
        const wall_rad = wall_rads[i];
        if (overlap_hole_region(wall_rad, coords1)) { 
            wall_list.push(false);
        }
        else if (overlap_hole_region(wall_rad, coords2)) {
            wall_list.push(false);
        }
        else {
            wall_list.push(true);
        }
    }
    return wall_list;
}

function circumference_coords(n_segs, r, ring_space, thickness, rot) {
    // get distance of start and end pos based on radius and ring_space
    // 1. calculate circumference
    var c = 2*Math.PI*r;
    // 2. standardize the ring_space
    var standard_ring_space = (ring_space + thickness)/c; // leave room for circular edges          thickness
    // 3. add correct measure for canvas arc method (2*pi)
    var end_start_dist = standard_ring_space*2*Math.PI;

    // get coordinates based on number of segments
    var step = 2*Math.PI/n_segs - end_start_dist; // - (n_segs*end_start_dist);
    // apply rotation (in percent of max)
    var rot_add = rot*2*Math.PI;
    var current_start = (end_start_dist/2 + rot_add)%(2*Math.PI); // so that first hole is centered
    var current_end = (current_start + step)%(2*Math.PI);
    var output_edges = [];          // for drawing the circular segments
    var output_ball_coords = [];    // for drawing the circular edges of the segments
    for (let index = 0; index < n_segs; index++) {
        // append to segment list
        output_edges.push({start: current_start, end: current_end});
        // append to edge coord list
        output_ball_coords.push(get_exact_coord(current_start, r, false));
        output_ball_coords.push(get_exact_coord(current_end, r, false));
        // step coords forward
        current_start = (current_end + end_start_dist)%(2*Math.PI);
        current_end = (current_start + step)%(2*Math.PI);
    }
    return [output_edges, output_ball_coords];
}

function circumference_coords_simpler(n_segs, rot, random=0) {
    // no pairs of coordinates, only one coordinate!
    // to get pairs, the circumference_coords_given_coords function needs to be called
    // apply rotation (in percent of max)

    var rot_add = rot*2*Math.PI;
    var current_start = rot_add%(2*Math.PI); // so that first hole is centered
    var output_edges = [];          // for drawing the circular segments

    for (let index = 0; index < n_segs; index++) {
        // append to segment list
        output_edges.push(current_start);
        current_start = (current_start + (2*Math.PI)/n_segs + random)%(2*Math.PI);
    }
    return output_edges;
}

function circumference_coords_rand(n_segs) {
    // apply rotation (in percent of max)
    var rot_add = rot*2*Math.PI;
    var current_start = rot_add; // so that first hole is centered
    var output_edges = [];          // for drawing the circular segments
    for (let index = 0; index < n_segs; index++) {
        // var rand_step = Math.random()*(rand_max - rand_min) + rand_min;
        // append to segment list
        output_edges.push(current_start);
        current_start = (current_start + (2*Math.PI))%(2*Math.PI);
    }
    return output_edges;
}

function circumference_coords_given_coords(coords, r, ring_space, thickness) {
    // get distance of start and end pos based on radius and ring_space
    // 1. calculate circumference
    var c = 2*Math.PI*r;
    // 2. standardize the ring_space
    var standard_ring_space = (ring_space + thickness)/c; // leave room for circular edges          thickness
    // 3. add correct measure for canvas arc method (2*pi)
    var end_start_dist = standard_ring_space*2*Math.PI;

    // get coordinates based on number of segments
    var step = 2*Math.PI/coords.length - end_start_dist; // - (n_segs*end_start_dist);
    var output_edges = [];          // for drawing the circular segments
    var output_ball_coords = [];    // for drawing the circular edges of the segments
    for (let index = 0; index < coords.length; index++) {
        var current_start = coords[index]%(2*Math.PI); // so that first hole is centered
        var current_end = (current_start + step)%(2*Math.PI);
        // append to segment list
        output_edges.push({start: current_start, end: current_end});
        // append to edge coord list
        output_ball_coords.push(get_exact_coord(current_start, r, false));
        output_ball_coords.push(get_exact_coord(current_end, r, false));
    }
    return [output_edges, output_ball_coords];
}

function circumference_from_diameter(d) {
    return Math.PI*d;
}

function start_end_from_coord(coord, c, ring_spacing, ring_thickness, flipped=false) {
    var standard_ring_space = (ring_spacing + ring_thickness)/c;
    var end_start_dist = standard_ring_space*2*Math.PI;
    // half spacing to left and to right
    var s = (coord - 0.5*end_start_dist + (2*Math.PI))%(2*Math.PI);
    var e = (coord + 0.5*end_start_dist + (2*Math.PI))%(2*Math.PI);
    if (flipped) {
        return {start: e, end: s};    
    } else {
        return {start: s, end: e};
    }
}

function equidistant_rads(n) {
    // centers first hole to top (i.e. 1.5*pi)
    var coords = [];
    var spacing = (2*Math.PI)/n;
    var current_pos = (1.5*Math.PI + ((Math.PI*2)/n/2))%(Math.PI*2); // 1.5*Math.PI;
    for (let index = 0; index < n; index++) {
        coords.push(current_pos);
        current_pos = (current_pos + spacing)%(2*Math.PI);
    }
    return coords
}

function equidistant_holes_coords(n, c, ring_spacing, ring_thickness) {
    // centers first hole to top (i.e. 1.5*pi)
    // 1. set all the relevant coords without start/end
    var coords = equidistant_rads(n);
    // 2. set the start end coords accordingly (acoordingly LOL)
    var output = [];
    for (let index = 0; index < n; index++) {
        output.push(start_end_from_coord(coords[index], c, ring_spacing, ring_thickness));
    }
    return output;
}

function equidistant_wall_coords(n) {
    var coords = [];
    var spacing = (2*Math.PI)/n;
    // set start position:
    var current_pos = 1.5*Math.PI; // (1.5*Math.PI + ((Math.PI*2)/n/2))%(Math.PI*2);
    // var current_pos = 1.5*Math.PI; // important issue: how to set first pos correctly
    for (let index = 0; index < n; index++) {
        coords.push(current_pos);
        current_pos = (current_pos + spacing)%(2*Math.PI);
    }
    return coords;
}