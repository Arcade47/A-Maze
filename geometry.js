var center_coord = {x: canvas.width/2, y: canvas.height/2};

function rad_to_degrees(rad) {
    var unitless = rad / 2*Math.PI;
    return unitless*360;
}

function degrees_to_rad(degree) {
    var unitless = degree / 360;
    return unitless*2*Math.PI;
}

/*
circ: {radius: ..., pos: {x: ..., y: ...}}
*/

function distance(pos1, pos2) {
    var xdiff = pos2.x - pos1.x;
    var ydiff = pos2.y - pos1.y;
    return Math.sqrt(xdiff*xdiff + ydiff*ydiff);
}

function overlap(circ1, circ2) {
    var sum_radii = circ1.radius + circ2.radius;
    if (distance(circ1.pos, circ2.pos) < sum_radii) {
        return true;
    } else {
        return false;
    }
}

function rad_to_coord(rad, distance) {
    var x_val = Math.cos(rad)*distance + center_coord.x;
    var y_val = Math.sin(rad)*distance + center_coord.y;
    return {x: x_val, y: y_val};
}

function coord_to_rad(coord) {
    var val = Math.atan2(coord.y - center_coord.y, coord.x - center_coord.x);
    if (val < 0) {
        return val + 2*Math.PI;
    } else {
        return val;
    }
}

function get_perimeter(radius) {
    return 2*Math.PI*radius;
}

function wrap_around(rad) {
    if (rad == 2*Math.PI) {
        return rad;
    }
    return rad%(2*Math.PI);
}

function round_down(number, size) {
    var n_fit = Math.floor(number/size);
    return n_fit*size;
}

function round_down_to_power_of_two(number) {
    var current_pot = 2;
    var n_fit = Math.floor(number/current_pot);
    while (n_fit > 0) {
        current_pot *= 2;
        n_fit = Math.floor(number/current_pot);
    }
    return current_pot/2;
}

function get_walls_rads(radius, thickness, spacing) {

    var output = [];

    var perimeter = get_perimeter(radius);

    // calculate number of walls
    // TODO: adjust
    var hole_wall_width = spacing + thickness;
    var max_walls = perimeter / hole_wall_width;
    max_walls = round_down_to_power_of_two(max_walls);

    // add the radians values
    var stepsize = (2*Math.PI)/max_walls;
    var current_rad = 1.5*Math.PI;
    for (let index = 0; index < max_walls; index++) {
        output.push(current_rad);
        current_rad += stepsize;
        current_rad = wrap_around(current_rad);
    }

    return output;
}

function rad_dist(radius, len) {
    var perimeter = get_perimeter(radius);
    return (len/perimeter)*(2*Math.PI);
}

function rad_between_rads(rad, rads) {
    if (rads[0] < rad && rad < rads[1]) {
        return true;
    } else if (rads[1] < rads[0]) {
        // end value is wrapped around maximum --> rads[1] < rads[0]
        var before_wrap = rads[0] < rad && rad <= (2*Math.PI);
        var after_wrap = rad < rads[1] && rad >= 0;
        if (before_wrap || after_wrap) {
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

function sort_list_to_another_list(l, other_l) {
    // zip together
    var zipper = l.map(function(e, i){ return [e, other_l[i]]; });
    // sort
    zipper.sort( function(a,b){ return a[1] < b[1]; });
    // remove second element ("other") of sublists
    return zipper.map(function(e, i){ return e[0];} );
}

function get_object_list_copy(l) {
    var output = [];
    for (let index = 0; index < l.length; index++) {
        const element = l[index];
        output.push(element);
    }
    return output;
}

// TODO: delete
function rad_from_ind(n_max_rads, ind_in_circ) {
    var first_center_rad = (1.5*Math.PI) + (2*Math.PI)/(n_max_rads*2);
    var stepsize = (2*Math.PI)/n_max_rads;
    return wrap_around(first_center_rad + (ind_in_circ*stepsize));
}
// function coord_from_rad_distance(coord, distance, rad) {
//     output_x = distance*Math.cos(rad);
//     output_y = distance*Math.sin(rad);
//     return {x: output_x + coord.x, y: output_y + coord.y};
// }