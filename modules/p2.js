
export class Line {
  /** Compute the intersection between two lines.
    * @param {Array}  l1          Line vector 1
   * @param {Array}  l2          Line vector 2
   * @param {Number} precision   Precision to use when checking if the lines are parallel
   * @returns {Array}              The intersection point.
   */
  static lineInt(l1, l2, precision) {
    precision = precision || 0;
    let i = [0, 0]; // point
    let a1, b1, c1, a2, b2, c2, det; // scalars
    a1 = l1[1][1] - l1[0][1];
    b1 = l1[0][0] - l1[1][0];
    c1 = a1 * l1[0][0] + b1 * l1[0][1];
    a2 = l2[1][1] - l2[0][1];
    b2 = l2[0][0] - l2[1][0];
    c2 = a2 * l2[0][0] + b2 * l2[0][1];
    det = a1 * b2 - a2 * b1;
    if (!Scalar.eq(det, 0, precision)) { // lines are not parallel
      i[0] = (b2 * c1 - b1 * c2) / det;
      i[1] = (a1 * c2 - a2 * c1) / det;
    }
    return i;
  }
  /**Checks if two line segments intersects.
   * @param {Array} p1 The start vertex of the first line segment.
   * @param {Array} p2 The end vertex of the first line segment.
   * @param {Array} q1 The start vertex of the second line segment.
   * @param {Array} q2 The end vertex of the second line segment.
   * @returns {Boolean} True if the two line segments intersect
   */
  static segmentsIntersect(p1, p2, q1, q2) {
    let dx = p2[0] - p1[0];
    let dy = p2[1] - p1[1];
    let da = q2[0] - q1[0];
    let db = q2[1] - q1[1];
    // segments are parallel
    if (da * dy - db * dx == 0)
      return false;
    let s = (dx * (q1[1] - p1[1]) + dy * (p1[0] - q1[0])) / (da * dy - db * dx)
    let t = (da * (p1[1] - q1[1]) + db * (q1[0] - p1[0])) / (db * dx - da * dy)
    return (s >= 0 && s <= 1 && t >= 0 && t <= 1);
  }
}

export class Point {
  /**Get the area of a triangle spanned by the three given points. Note that the area will be negative if the points are not given in counter-clockwise order.
    * @param {Array} a
   * @param {Array} b
   * @param {Array} c
   * @returns {Number}
   */
  static area(a, b, c) {
    return (((b[0] - a[0]) * (c[1] - a[1])) - ((c[0] - a[0]) * (b[1] - a[1])));
  }
  left(a, b, c) {
    return Point.area(a, b, c) > 0;
  }
  leftOn(a, b, c) {
    return Point.area(a, b, c) >= 0;
  }
  right(a, b, c) {
    return Point.area(a, b, c) < 0;
  }
  rightOn(a, b, c) {
    return Point.area(a, b, c) <= 0;
  }
  static tmpPoint1 = []
  static tmpPoint2 = [];
  /**Check if three points are collinear
   * @param {Array} a
   * @param {Array} b
   * @param {Array} c
   * @param {Number} [thresholdAngle=0] Threshold angle to use when comparing the vectors. The function will return true if the angle between the resulting vectors is less than this value. Use zero for max precision.
   * @returns {Boolean}
   */
  collinear(a, b, c, thresholdAngle) {
    if (!thresholdAngle)
      return Point.area(a, b, c) == 0;
    else {
      let ab = Point.tmpPoint1,
        bc = Point.tmpPoint2;
      ab[0] = b[0] - a[0];
      ab[1] = b[1] - a[1];
      bc[0] = c[0] - b[0];
      bc[1] = c[1] - b[1];
      let dot = ab[0] * bc[0] + ab[1] * bc[1],
        magA = Math.sqrt(ab[0] * ab[0] + ab[1] * ab[1]),
        magB = Math.sqrt(bc[0] * bc[0] + bc[1] * bc[1]),
        angle = Math.acos(dot / (magA * magB));
      return angle < thresholdAngle;
    }
  }
  sqdist(a, b) {
    let dx = b[0] - a[0];
    let dy = b[1] - a[1];
    return dx * dx + dy * dy;
  }
}

export class Polygon {
  constructor() {
    /**@type {Array} */
    this.vertices = new Array();
  }
  /**Get a vertex at position i. It does not matter if i is out of bounds, this function will just cycle.
   * @param {Number} i
   * @returns {Array}
   */
  at(i) {
    let v = this.vertices,
      s = v.length;
    return v[i < 0 ? i % s + s : i % s];
  }
  /**Get first vertex
   * @returns {Array}
   */
  first() {
    return this.vertices[0];
  }
  /**Get last vertex
   * @returns {Array}
   */
  last() {
    return this.vertices[this.vertices.length - 1];
  }
  /**Clear the polygon data*/
  clear() {
    this.vertices.length = 0;
  }
  /**Append points "from" to "to"-1 from an other polygon "poly" onto this one.
   * @param {Polygon} poly The polygon to get points from.
   * @param {Number} from The vertex index in "poly".
   * @param {Number} to The end vertex index in "poly". Note that this vertex is NOT included when appending.
   */
  append(poly, from = undefined, to = undefined) {
    if (from === undefined) from = 0;
    if (to === undefined) to = poly.vertices.length;
    for (let i = from; i < to; i++) {
      this.vertices.push(poly.vertices[i]);
    }
  }
  /**Make sure that the polygon vertices are ordered counter-clockwise.
   */
  makeCCW() {
    let br = 0,
      v = this.vertices;
    // find bottom right point
    for (let i = 1; i < this.vertices.length; ++i) {
      if (v[i][1] < v[br][1] || (v[i][1] == v[br][1] && v[i][0] > v[br][0])) {
        br = i;
      }
    }
    // reverse poly if clockwise
    if (!Point.left(this.at(br - 1), this.at(br), this.at(br + 1))) {
      this.reverse();
    }
  }
  /**Reverse the vertices in the polygon
   */
  reverse() {
    this.vertices.reverse();
  }
  /**Check if a point in the polygon is a reflex point
   * @param {Number}  i
   * @returns {Boolean}
   */
  isReflex(i) {
    return Point.right(this.at(i - 1), this.at(i), this.at(i + 1));
  }
  static tmpLine1 = []
  static tmpLine2 = [];
  /**Check if two vertices in the polygon can see each other
   * @param {Number} a Vertex index 1
   * @param {Number} b Vertex index 2
   * @returns {Boolean}
   */
  canSee(a, b) {
    let p, dist, l1 = tmpLine1, l2 = tmpLine2;
    if (Point.leftOn(this.at(a + 1), this.at(a), this.at(b)) && Point.rightOn(this.at(a - 1), this.at(a), this.at(b))) {
      return false;
    }
    dist = Point.sqdist(this.at(a), this.at(b));
    for (let i = 0; i !== this.vertices.length; ++i) { // for each edge
      if ((i + 1) % this.vertices.length === a || i === a) // ignore incident edges
        continue;
      if (Point.leftOn(this.at(a), this.at(b), this.at(i + 1)) && Point.rightOn(this.at(a), this.at(b), this.at(i))) { // if diag intersects an edge
        l1[0] = this.at(a);
        l1[1] = this.at(b);
        l2[0] = this.at(i);
        l2[1] = this.at(i + 1);
        p = Line.lineInt(l1, l2);
        if (Point.sqdist(this.at(a), p) < dist) { // if edge is blocking visibility to b
          return false;
        }
      }
    }
    return true;
  }
  /** Copy the polygon from vertex i to vertex j.
   * @param {Number} i
   * @param {Number} j
   * @param {Polygon} [targetPoly]   Optional target polygon to save in.
   * @returns {Polygon} The resulting copy.
   */
  copy(i, j, targetPoly) {
    let p = targetPoly || new Polygon();
    p.clear();
    if (i < j) {
      // Insert all vertices from i to j
      for (let k = i; k <= j; k++)
        p.vertices.push(this.vertices[k]);
    } else {
      // Insert vertices 0 to j
      for (let k = 0; k <= j; k++)
        p.vertices.push(this.vertices[k]);
      // Insert vertices i to end
      for (let k = i; k < this.vertices.length; k++)
        p.vertices.push(this.vertices[k]);
    }
    return p;
  }
  /**Decomposes the polygon into convex pieces. Returns a list of edges [[p1,p2],[p2,p3],...] that cuts the polygon.
   * Note that this algorithm has complexity O(N^4) and will be very slow for polygons with many vertices.
   * @returns {Array}
   */
  getCutEdges() {
    let min = [], tmp1 = [], tmp2 = [], tmpPoly = new Polygon();
    let nDiags = Number.MAX_VALUE;
    for (let i = 0; i < this.vertices.length; ++i) {
      if (this.isReflex(i)) {
        for (let j = 0; j < this.vertices.length; ++j) {
          if (this.canSee(i, j)) {
            tmp1 = this.copy(i, j, tmpPoly).getCutEdges();
            tmp2 = this.copy(j, i, tmpPoly).getCutEdges();
            for (let k = 0; k < tmp2.length; k++)
              tmp1.push(tmp2[k]);
            if (tmp1.length < nDiags) {
              min = tmp1;
              nDiags = tmp1.length;
              min.push([this.at(i), this.at(j)]);
            }
          }
        }
      }
    }
    return min;
  }
  /**Decomposes the polygon into one or more convex sub-Polygons.
   * @returns {Array} An array or Polygon objects.
   */
  decomp() {
    let edges = this.getCutEdges();
    if (edges.length > 0)
      return this.slice(edges);
    else
      return [this];
  }
  /**Slices the polygon given one or more cut edges. If given one, this function will return two polygons (false on failure). If many, an array of polygons.
   * @param {Array} cutEdges A list of edges, as returned by .getCutEdges()
   * @returns {Array}
   */
  slice(cutEdges) {
    if (cutEdges.length == 0) return [this];
    if (cutEdges instanceof Array && cutEdges.length && cutEdges[0] instanceof Array && cutEdges[0].length == 2 && cutEdges[0][0] instanceof Array) {
      let polys = [this];
      for (let i = 0; i < cutEdges.length; i++) {
        let cutEdge = cutEdges[i];
        // Cut all polys
        for (let j = 0; j < polys.length; j++) {
          let poly = polys[j];
          let result = poly.slice(cutEdge);
          if (result) {
            // Found poly! Cut and quit
            polys.splice(j, 1);
            polys.push(result[0], result[1]);
            break;
          }
        }
      }
      return polys;
    } else {
      // Was given one edge
      let cutEdge = cutEdges;
      let i = this.vertices.indexOf(cutEdge[0]);
      let j = this.vertices.indexOf(cutEdge[1]);
      if (i != -1 && j != -1) {
        return [this.copy(i, j),
        this.copy(j, i)];
      } else {
        return false;
      }
    }
  }
  /**Checks that the line segments of this polygon do not intersect each other.
   * @param {Array} path An array of vertices e.g. [[0,0],[0,1],...]
   * @returns {Boolean}
   * @todo Should it check all segments with all others?
   */
  isSimple() {
    let path = this.vertices;
    // Check
    for (let i = 0; i < path.length - 1; i++) {
      for (let j = 0; j < i - 1; j++) {
        if (Line.segmentsIntersect(path[i], path[i + 1], path[j], path[j + 1])) {
          return false;
        }
      }
    }
    // Check the segment between the last and the first point to all others
    for (let i = 1; i < path.length - 2; i++) {
      if (Line.segmentsIntersect(path[0], path[path.length - 1], path[i], path[i + 1])) {
        return false;
      }
    }
    return true;
  }
  static getIntersectionPoint(p1, p2, q1, q2, delta) {
    delta = delta || 0;
    let a1 = p2[1] - p1[1];
    let b1 = p1[0] - p2[0];
    let c1 = (a1 * p1[0]) + (b1 * p1[1]);
    let a2 = q2[1] - q1[1];
    let b2 = q1[0] - q2[0];
    let c2 = (a2 * q1[0]) + (b2 * q1[1]);
    let det = (a1 * b2) - (a2 * b1);
    if (!Scalar.eq(det, 0, delta))
      return [((b2 * c1) - (b1 * c2)) / det, ((a1 * c2) - (a2 * c1)) / det]
    else
      return [0, 0]
  }
  /**Quickly decompose the Polygon into convex sub-polygons.
   * @param {Array} result
   * @param {Array} reflexVertices
   * @param {Array} steinerPoints
   * @param {Number} delta
   * @param {Number} maxlevel
   * @param {Number} level
   * @returns {Array}
   */
  quickDecomp(result, reflexVertices, steinerPoints, delta = 25, maxlevel = 100, level = 0) {
    result = typeof (result) != "undefined" ? result : [];
    reflexVertices = reflexVertices || [];
    steinerPoints = steinerPoints || [];
    let upperInt = [0, 0], lowerInt = [0, 0], p = [0, 0]; // Points
    let upperDist = 0, lowerDist = 0, d = 0, closestDist = 0; // scalars
    let upperIndex = 0, lowerIndex = 0, closestIndex = 0; // Integers
    let lowerPoly = new Polygon(), upperPoly = new Polygon(); // polygons
    let poly = this,
      v = this.vertices;
    if (v.length < 3) return result;
    level++;
    if (level > maxlevel) {
      console.warn("quickDecomp: max level (" + maxlevel + ") reached.");
      return result;
    }
    for (let i = 0; i < this.vertices.length; ++i) {
      if (poly.isReflex(i)) {
        reflexVertices.push(poly.vertices[i]);
        upperDist = lowerDist = Number.MAX_VALUE;

        for (let j = 0; j < this.vertices.length; ++j) {
          if (Point.left(poly.at(i - 1), poly.at(i), poly.at(j))
            && Point.rightOn(poly.at(i - 1), poly.at(i), poly.at(j - 1))) { // if line intersects with an edge
            p = getIntersectionPoint(poly.at(i - 1), poly.at(i), poly.at(j), poly.at(j - 1)); // find the point of intersection
            if (Point.right(poly.at(i + 1), poly.at(i), p)) { // make sure it's inside the poly
              d = Point.sqdist(poly.vertices[i], p);
              if (d < lowerDist) { // keep only the closest intersection
                lowerDist = d;
                lowerInt = p;
                lowerIndex = j;
              }
            }
          }
          if (Point.left(poly.at(i + 1), poly.at(i), poly.at(j + 1))
            && Point.rightOn(poly.at(i + 1), poly.at(i), poly.at(j))) {
            p = getIntersectionPoint(poly.at(i + 1), poly.at(i), poly.at(j), poly.at(j + 1));
            if (Point.left(poly.at(i - 1), poly.at(i), p)) {
              d = Point.sqdist(poly.vertices[i], p);
              if (d < upperDist) {
                upperDist = d;
                upperInt = p;
                upperIndex = j;
              }
            }
          }
        }
        // if there are no vertices to connect to, choose a point in the middle
        if (lowerIndex == (upperIndex + 1) % this.vertices.length) {
          //console.log("Case 1: Vertex("+i+"), lowerIndex("+lowerIndex+"), upperIndex("+upperIndex+"), poly.size("+this.vertices.length+")");
          p[0] = (lowerInt[0] + upperInt[0]) / 2;
          p[1] = (lowerInt[1] + upperInt[1]) / 2;
          steinerPoints.push(p);
          if (i < upperIndex) {
            //lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.begin() + upperIndex + 1);
            lowerPoly.append(poly, i, upperIndex + 1);
            lowerPoly.vertices.push(p);
            upperPoly.vertices.push(p);
            if (lowerIndex != 0) {
              //upperPoly.insert(upperPoly.end(), poly.begin() + lowerIndex, poly.end());
              upperPoly.append(poly, lowerIndex, poly.vertices.length);
            }
            //upperPoly.insert(upperPoly.end(), poly.begin(), poly.begin() + i + 1);
            upperPoly.append(poly, 0, i + 1);
          } else {
            if (i != 0) {
              //lowerPoly.insert(lowerPoly.end(), poly.begin() + i, poly.end());
              lowerPoly.append(poly, i, poly.vertices.length);
            }
            //lowerPoly.insert(lowerPoly.end(), poly.begin(), poly.begin() + upperIndex + 1);
            lowerPoly.append(poly, 0, upperIndex + 1);
            lowerPoly.vertices.push(p);
            upperPoly.vertices.push(p);
            //upperPoly.insert(upperPoly.end(), poly.begin() + lowerIndex, poly.begin() + i + 1);
            upperPoly.append(poly, lowerIndex, i + 1);
          }
        } else {
          // connect to the closest point within the triangle
          //console.log("Case 2: Vertex("+i+"), closestIndex("+closestIndex+"), poly.size("+this.vertices.length+")\n");
          if (lowerIndex > upperIndex) {
            upperIndex += this.vertices.length;
          }
          closestDist = Number.MAX_VALUE;
          if (upperIndex < lowerIndex) {
            return result;
          }
          for (let j = lowerIndex; j <= upperIndex; ++j) {
            if (Point.leftOn(poly.at(i - 1), poly.at(i), poly.at(j))
              && Point.rightOn(poly.at(i + 1), poly.at(i), poly.at(j))) {
              d = Point.sqdist(poly.at(i), poly.at(j));
              if (d < closestDist) {
                closestDist = d;
                closestIndex = j % this.vertices.length;
              }
            }
          }
          if (i < closestIndex) {
            lowerPoly.append(poly, i, closestIndex + 1);
            if (closestIndex != 0) {
              upperPoly.append(poly, closestIndex, v.length);
            }
            upperPoly.append(poly, 0, i + 1);
          } else {
            if (i != 0) {
              lowerPoly.append(poly, i, v.length);
            }
            lowerPoly.append(poly, 0, closestIndex + 1);
            upperPoly.append(poly, closestIndex, i + 1);
          }
        }
        // solve smallest poly first
        if (lowerPoly.vertices.length < upperPoly.vertices.length) {
          lowerPoly.quickDecomp(result, reflexVertices, steinerPoints, delta, maxlevel, level);
          upperPoly.quickDecomp(result, reflexVertices, steinerPoints, delta, maxlevel, level);
        } else {
          upperPoly.quickDecomp(result, reflexVertices, steinerPoints, delta, maxlevel, level);
          lowerPoly.quickDecomp(result, reflexVertices, steinerPoints, delta, maxlevel, level);
        }
        return result;
      }
    }
    result.push(this);
    return result;
  }
  /**Remove collinear points in the polygon.
   * @param {Number} precision The threshold angle to use when determining whether two edges are collinear. Use zero for finest precision.
   * @returns {Number}           The number of points removed
   */
  removeCollinearPoints(precision) {
    let num = 0;
    for (let i = this.vertices.length - 1; this.vertices.length > 3 && i >= 0; --i) {
      if (Point.collinear(this.at(i - 1), this.at(i), this.at(i + 1), precision)) {
        // Remove the middle point
        this.vertices.splice(i % this.vertices.length, 1);
        i--; // Jump one point forward. Otherwise we may get a chain removal
        num++;
      }
    }
    return num;
  }
}
export class Scalar {
  /**Check if two scalars are equal
   * @param {Number} a
   * @param {Number} b
   * @param {Number} [precision]
   * @returns {Boolean}
   */
  static eq(a, b, precision = 0) {
    return Math.abs(a - b) < precision;
  }
}

export class AABB {
  /**Axis aligned bounding box class.
   * @param {{upperBound:Array, lowerBound:Array}} options
   */
  constructor(options) {
    /**@type {Array}*/
    this.lowerBound = vec2.create();
    if (options && options.lowerBound) {
      vec2.copy(this.lowerBound, options.lowerBound);
    }
    /**The upper bound of the bounding box.
    * @property upperBound
    * @type {Array}
    */
    this.upperBound = vec2.create();
    if (options && options.upperBound) {
      vec2.copy(this.upperBound, options.upperBound);
    }
  }
  static tmp = vec2.create();
  /**Set the AABB bounds from a set of points, transformed by the given position and angle.
   * @param {Array} points An array of vec2's.
   * @param {Array} position
   * @param {number} angle
   * @param {number} skinSize Some margin to be added to the AABB.
   */
  setFromPoints(points, position, angle, skinSize) {
    let l = this.lowerBound,
      u = this.upperBound;
    if (typeof (angle) !== "number") {
      angle = 0;
    }
    // Set to the first point
    if (angle !== 0) {
      vec2.rotate(l, points[0], angle);
    } else {
      vec2.copy(l, points[0]);
    }
    vec2.copy(u, l);
    // Compute cosines and sines just once
    let cosAngle = Math.cos(angle),
      sinAngle = Math.sin(angle);
    for (let i = 1; i < points.length; i++) {
      let p = points[i];
      if (angle !== 0) {
        let x = p[0],
          y = p[1];
        AABB.tmp[0] = cosAngle * x - sinAngle * y;
        AABB.tmp[1] = sinAngle * x + cosAngle * y;
        p = AABB.tmp;
      }
      for (let j = 0; j < 2; j++) {
        if (p[j] > u[j]) {
          u[j] = p[j];
        }
        if (p[j] < l[j]) {
          l[j] = p[j];
        }
      }
    }
    // Add offset
    if (position) {
      vec2.add(this.lowerBound, this.lowerBound, position);
      vec2.add(this.upperBound, this.upperBound, position);
    }
    if (skinSize) {
      this.lowerBound[0] -= skinSize;
      this.lowerBound[1] -= skinSize;
      this.upperBound[0] += skinSize;
      this.upperBound[1] += skinSize;
    }
  }
  /**Copy bounds from an AABB to this AABB
   * @param {AABB} aabb
   */
  copy(aabb) {
    vec2.copy(this.lowerBound, aabb.lowerBound);
    vec2.copy(this.upperBound, aabb.upperBound);
  }
  /**Extend this AABB so that it covers the given AABB too.
   * @param {AABB} aabb
   */
  extend(aabb) {
    // Loop over x and y
    let i = 2;
    while (i--) {
      // Extend lower bound
      let l = aabb.lowerBound[i];
      if (this.lowerBound[i] > l) {
        this.lowerBound[i] = l;
      }
      // Upper
      let u = aabb.upperBound[i];
      if (this.upperBound[i] < u) {
        this.upperBound[i] = u;
      }
    }
  }
  /**Returns true if the given AABB overlaps this AABB.
   * @param {AABB} aabb
   * @returns {Boolean}
   */
  overlaps(aabb) {
    let l1 = this.lowerBound,
      u1 = this.upperBound,
      l2 = aabb.lowerBound,
      u2 = aabb.upperBound;
    //      l2        u2
    //      |---------|
    // |--------|
    // l1       u1
    return ((l2[0] <= u1[0] && u1[0] <= u2[0]) || (l1[0] <= u2[0] && u2[0] <= u1[0])) &&
      ((l2[1] <= u1[1] && u1[1] <= u2[1]) || (l1[1] <= u2[1] && u2[1] <= u1[1]));
  }
  /**@param {Array} point
   * @returns {Boolean}
   */
  containsPoint(point) {
    let l = this.lowerBound,
      u = this.upperBound;
    return l[0] <= point[0] && point[0] <= u[0] && l[1] <= point[1] && point[1] <= u[1];
  }
  /**Check if the AABB is hit by a ray.
   * @param {Ray} ray
   * @returns {number} -1 if no hit, a number between 0 and 1 if hit.
   */
  overlapsRay(ray) {
    let t = 0;
    // ray.direction is unit direction vector of ray
    let dirFracX = 1 / ray.direction[0];
    let dirFracY = 1 / ray.direction[1];
    // this.lowerBound is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
    let t1 = (this.lowerBound[0] - ray.from[0]) * dirFracX;
    let t2 = (this.upperBound[0] - ray.from[0]) * dirFracX;
    let t3 = (this.lowerBound[1] - ray.from[1]) * dirFracY;
    let t4 = (this.upperBound[1] - ray.from[1]) * dirFracY;
    let tmin = Math.max(Math.max(Math.min(t1, t2), Math.min(t3, t4)));
    let tmax = Math.min(Math.min(Math.max(t1, t2), Math.max(t3, t4)));
    // if tmax < 0, ray (line) is intersecting AABB, but whole AABB is behing us
    if (tmax < 0) {
      //t = tmax;
      return -1;
    }
    // if tmin > tmax, ray doesn't intersect AABB
    if (tmin > tmax) {
      //t = tmax;
      return -1;
    }
    return tmin;
  }
}

export class Broadphase {
  constructor(type) {
    this.type = type;
    /**@type {Array}*/
    this.result = new Array();
    /**@type {World}*/
    this.world = null;
    this.boundingVolumeType = Broadphase.AABB;
  }
  /**Axis aligned bounding box type.
   */
  static AABB = 1;
  static BOUNDING_CIRCLE = 2;
  /**Set the world that we are searching for collision pairs in
   * @param {World} world
   */
  setWorld(world) {
    this.world = world;
  }
  /**Get all potential intersecting body pairs.
   * @param {World} world The world to search in.
   * @returns {Array} An array of the bodies, ordered in pairs. Example: A result of [a,b,c,d] means that the potential pairs are: (a,b), (c,d).
   */
  getCollisionPairs(world) { }
  static dist = vec2.create();
  /**Check whether the bounding radius of two bodies overlap.
   * @param {Body} bodyA
   * @param {Body} bodyB
   * @returns {Boolean}
   */
  boundingRadiusCheck(bodyA, bodyB) {
    vec2.sub(Broadphase.dist, bodyA.position, bodyB.position);
    let d2 = vec2.squaredLength(Broadphase.dist),
      r = bodyA.boundingRadius + bodyB.boundingRadius;
    return d2 <= r * r;
  }
  /**Check whether the bounding radius of two bodies overlap.
   * @param {Body} bodyA
   * @param {Body} bodyB
   * @returns {Boolean}
   */
  aabbCheck(bodyA, bodyB) {
    return bodyA.getAABB().overlaps(bodyB.getAABB());
  }
  /**Check whether the bounding radius of two bodies overlap.
   * @param {Body} bodyA
   * @param {Body} bodyB
   * @returns {Boolean}
   */
  boundingVolumeCheck(bodyA, bodyB) {
    let result;
    switch (this.boundingVolumeType) {
      case Broadphase.BOUNDING_CIRCLE:
        result = Broadphase.boundingRadiusCheck(bodyA, bodyB);
        break;
      case Broadphase.AABB:
        result = Broadphase.aabbCheck(bodyA, bodyB);
        break;
      default:
        throw new Error('Bounding volume type not recognized: ' + this.boundingVolumeType);
    }
    return result;
  }
  /**Check whether two bodies are allowed to collide at all.
   * @param {Body} bodyA
   * @param {Body} bodyB
   * @returns {Boolean}
   */
  canCollide(bodyA, bodyB) {
    let KINEMATIC = Body.KINEMATIC;
    let STATIC = Body.STATIC;
    // Cannot collide static bodies
    if (bodyA.type === STATIC && bodyB.type === STATIC) {
      return false;
    }
    // Cannot collide static vs kinematic bodies
    if ((bodyA.type === KINEMATIC && bodyB.type === STATIC) ||
      (bodyA.type === STATIC && bodyB.type === KINEMATIC)) {
      return false;
    }
    // Cannot collide kinematic vs kinematic
    if (bodyA.type === KINEMATIC && bodyB.type === KINEMATIC) {
      return false;
    }
    // Cannot collide both sleeping bodies
    if (bodyA.sleepState === Body.SLEEPING && bodyB.sleepState === Body.SLEEPING) {
      return false;
    }
    // Cannot collide if one is static and the other is sleeping
    if ((bodyA.sleepState === Body.SLEEPING && bodyB.type === STATIC) ||
      (bodyB.sleepState === Body.SLEEPING && bodyA.type === STATIC)) {
      return false;
    }
    return true;
  }
  static NAIVE = 1;
  static SAP = 2;
}

/**Naive broadphase implementation. Does N^2 tests*/
export class NaiveBroadphase extends Broadphase {
  constructor() {
    super(Broadphase.NAIVE);
  }
  /**Get the colliding pairs
   * @param {World} world
   * @returns {Array}
   */
  getCollisionPairs(world) {
    let bodies = world.bodies,
      result = this.result;
    result.length = 0;
    for (let i = 0, Ncolliding = bodies.length; i !== Ncolliding; i++) {
      let bi = bodies[i];
      for (let j = 0; j < i; j++) {
        let bj = bodies[j];
        if (Broadphase.canCollide(bi, bj) && this.boundingVolumeCheck(bi, bj)) {
          result.push(bi, bj);
        }
      }
    }
    return result;
  }
  /**Returns all the bodies within an AABB.
   * @param {World} world
   * @param {AABB} aabb
   * @param {array} result An array to store resulting bodies in.
   * @returns {array}
   */
  aabbQuery(world, aabb, result) {
    result = result || [];
    let bodies = world.bodies;
    for (let i = 0; i < bodies.length; i++) {
      let b = bodies[i];
      if (b.aabbNeedsUpdate) {
        b.updateAABB();
      }
      if (b.aabb.overlaps(aabb)) {
        result.push(b);
      }
    }
    return result;
  }
}
let yAxis = vec2.fromValues(0, 1);
let tmp1 = vec2.fromValues(0, 0)
  , tmp2 = vec2.fromValues(0, 0)
  , tmp3 = vec2.fromValues(0, 0)
  , tmp4 = vec2.fromValues(0, 0)
  , tmp5 = vec2.fromValues(0, 0)
  , tmp6 = vec2.fromValues(0, 0)
  , tmp7 = vec2.fromValues(0, 0)
  , tmp8 = vec2.fromValues(0, 0)
  , tmp9 = vec2.fromValues(0, 0)
  , tmp10 = vec2.fromValues(0, 0)
  , tmp11 = vec2.fromValues(0, 0)
  , tmp12 = vec2.fromValues(0, 0)
  , tmp13 = vec2.fromValues(0, 0)
  , tmp14 = vec2.fromValues(0, 0)
  , tmp15 = vec2.fromValues(0, 0)
  , tmp16 = vec2.fromValues(0, 0)
  , tmp17 = vec2.fromValues(0, 0)
  , tmp18 = vec2.fromValues(0, 0)
  , tmpArray = [];
/**Narrowphase. Creates contacts and friction given shapes and transforms
 */
export class Narrowphase {
  constructor() {
    this.contactEquations = new Array();
    this.frictionEquations = new Array();
    this.enableFriction = true;
    this.enabledEquations = true;
    this.slipForce = 10.0;
    this.frictionCoefficient = 0.3;
    this.surfaceVelocity = 0;
    this.contactEquationPool = new ContactEquationPool({ size: 32 });
    this.frictionEquationPool = new FrictionEquationPool({ size: 64 });
    this.restitution = 0;
    this.stiffness = Equation.DEFAULT_STIFFNESS;
    this.relaxation = Equation.DEFAULT_RELAXATION;
    this.frictionStiffness = Equation.DEFAULT_STIFFNESS;
    this.frictionRelaxation = Equation.DEFAULT_RELAXATION;
    this.enableFrictionReduction = true;
    this.collidingBodiesLastStep = new TupleDictionary();
    this.contactSkinSize = 0.01;
  }
  static bodiesOverlap_shapePositionA = vec2.create();
  static bodiesOverlap_shapePositionB = vec2.create();
  /**@param {Body} bodyA
   * @param {Body} bodyB
   * @returns {Boolean}
   * @todo shape world transforms are wrong
   */
  bodiesOverlap(bodyA, bodyB) {
    let shapePositionA = Narrowphase.bodiesOverlap_shapePositionA;
    let shapePositionB = Narrowphase.bodiesOverlap_shapePositionB;
    // Loop over all shapes of bodyA
    for (let k = 0, Nshapesi = bodyA.shapes.length; k !== Nshapesi; k++) {
      let shapeA = bodyA.shapes[k];
      bodyA.toWorldFrame(shapePositionA, shapeA.position);
      // All shapes of body j
      for (let l = 0, Nshapesj = bodyB.shapes.length; l !== Nshapesj; l++) {
        let shapeB = bodyB.shapes[l];
        bodyB.toWorldFrame(shapePositionB, shapeB.position);
        if (this[shapeA.type | shapeB.type](
          bodyA,
          shapeA,
          shapePositionA,
          shapeA.angle + bodyA.angle,
          bodyB,
          shapeB,
          shapePositionB,
          shapeB.angle + bodyB.angle,
          true
        )) {
          return true;
        }
      }
    }
    return false;
  }
  /**Check if the bodies were in contact since the last reset().
   * @param {Body} bodyA
   * @param {Body} bodyB
   * @returns {Boolean}
   */
  collidedLastStep(bodyA, bodyB) {
    let id1 = bodyA.id | 0,
      id2 = bodyB.id | 0;
    return this.collidingBodiesLastStep.get(id1, id2);
  }
  /**Throws away the old equations and gets ready to create new*/
  reset() {
    this.collidingBodiesLastStep.reset();
    let eqs = this.contactEquations;
    let l = eqs.length;
    while (l--) {
      let eq = eqs[l],
        id1 = eq.bodyA.id,
        id2 = eq.bodyB.id;
      this.collidingBodiesLastStep.set(id1, id2, true);
    }
    let ce = this.contactEquations,
      fe = this.frictionEquations;
    for (let i = 0; i < ce.length; i++) {
      this.contactEquationPool.release(ce[i]);
    }
    for (let i = 0; i < fe.length; i++) {
      this.frictionEquationPool.release(fe[i]);
    }
    // Reset
    this.contactEquations.length = this.frictionEquations.length = 0;
  }
  /**Creates a ContactEquation, either by reusing an existing object or creating a new one.
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @returns {ContactEquation}
  */
  createContactEquation(bodyA, bodyB, shapeA, shapeB) {
    let c = this.contactEquationPool.get();
    c.bodyA = bodyA;
    c.bodyB = bodyB;
    c.shapeA = shapeA;
    c.shapeB = shapeB;
    c.restitution = this.restitution;
    c.firstImpact = !this.collidedLastStep(bodyA, bodyB);
    c.stiffness = this.stiffness;
    c.relaxation = this.relaxation;
    c.needsUpdate = true;
    c.enabled = this.enabledEquations;
    c.offset = this.contactSkinSize;
    return c;
  }
  /**Creates a FrictionEquation, either by reusing an existing object or creating a new one.
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @returns {FrictionEquation}
  */
  createFrictionEquation(bodyA, bodyB, shapeA, shapeB) {
    let c = this.frictionEquationPool.get();
    c.bodyA = bodyA;
    c.bodyB = bodyB;
    c.shapeA = shapeA;
    c.shapeB = shapeB;
    c.setSlipForce(this.slipForce);
    c.frictionCoefficient = this.frictionCoefficient;
    c.relativeVelocity = this.surfaceVelocity;
    c.enabled = this.enabledEquations;
    c.needsUpdate = true;
    c.stiffness = this.frictionStiffness;
    c.relaxation = this.frictionRelaxation;
    c.contactEquations.length = 0;
    return c;
  }
  /**Creates a FrictionEquation given the data in the ContactEquation. Uses same offset vectors ri and rj, but the tangent vector will be constructed from the collision normal.
  * @param {ContactEquation} contactEquation
  * @returns {FrictionEquation}
  */
  createFrictionFromContact(c) {
    let eq = this.createFrictionEquation(c.bodyA, c.bodyB, c.shapeA, c.shapeB);
    vec2.copy(eq.contactPointA, c.contactPointA);
    vec2.copy(eq.contactPointB, c.contactPointB);
    vec2.rotate90cw(eq.t, c.normalA);
    eq.contactEquations.push(c);
    return eq;
  }
  // Take the average N latest contact point on the plane.
  createFrictionFromAverage(numContacts) {
    let c = this.contactEquations[this.contactEquations.length - 1];
    let eq = this.createFrictionEquation(c.bodyA, c.bodyB, c.shapeA, c.shapeB);
    let bodyA = c.bodyA;
    let bodyB = c.bodyB;
    vec2.set(eq.contactPointA, 0, 0);
    vec2.set(eq.contactPointB, 0, 0);
    vec2.set(eq.t, 0, 0);
    for (let i = 0; i !== numContacts; i++) {
      c = this.contactEquations[this.contactEquations.length - 1 - i];
      if (c.bodyA === bodyA) {
        vec2.add(eq.t, eq.t, c.normalA);
        vec2.add(eq.contactPointA, eq.contactPointA, c.contactPointA);
        vec2.add(eq.contactPointB, eq.contactPointB, c.contactPointB);
      } else {
        vec2.sub(eq.t, eq.t, c.normalA);
        vec2.add(eq.contactPointA, eq.contactPointA, c.contactPointB);
        vec2.add(eq.contactPointB, eq.contactPointB, c.contactPointA);
      }
      eq.contactEquations.push(c);
    }
    let invNumContacts = 1 / numContacts;
    vec2.scale(eq.contactPointA, eq.contactPointA, invNumContacts);
    vec2.scale(eq.contactPointB, eq.contactPointB, invNumContacts);
    vec2.normalize(eq.t, eq.t);
    vec2.rotate90cw(eq.t, eq.t);
    return eq;
  }
  /**Convex/line narrowphase
   * @param  {Body}       convexBody
   * @param  {Convex}     convexShape
   * @param  {Array}      convexOffset
   * @param  {Number}     convexAngle
   * @param  {Body}       lineBody
   * @param  {Line}       lineShape
   * @param  {Array}      lineOffset
   * @param  {Number}     lineAngle
   * @param {boolean}     justTest
   * @todo Implement me!
   */
  convexLine(
    convexBody,
    convexShape,
    convexOffset,
    convexAngle,
    lineBody,
    lineShape,
    lineOffset,
    lineAngle,
    justTest
  ) {
    // TODO
    if (justTest) {
      return false;
    } else {
      return 0;
    }
  }
  /**
 * Line/box narrowphase
 * @method lineBox
 * @param  {Body}       lineBody
 * @param  {Line}       lineShape
 * @param  {Array}      lineOffset
 * @param  {Number}     lineAngle
 * @param  {Body}       boxBody
 * @param  {Box}  boxShape
 * @param  {Array}      boxOffset
 * @param  {Number}     boxAngle
 * @param  {Boolean}    justTest
 * @todo Implement me!
 */
  lineBox(
    lineBody,
    lineShape,
    lineOffset,
    lineAngle,
    boxBody,
    boxShape,
    boxOffset,
    boxAngle,
    justTest
  ) {
    // TODO
    if (justTest) {
      return false;
    } else {
      return 0;
    }
  }

  static setConvexToCapsuleShapeMiddle(convexShape, capsuleShape) {
    vec2.set(convexShape.vertices[0], -capsuleShape.length * 0.5, -capsuleShape.radius);
    vec2.set(convexShape.vertices[1], capsuleShape.length * 0.5, -capsuleShape.radius);
    vec2.set(convexShape.vertices[2], capsuleShape.length * 0.5, capsuleShape.radius);
    vec2.set(convexShape.vertices[3], -capsuleShape.length * 0.5, capsuleShape.radius);
  }
  static convexCapsule_tempRect = new Box({ width: 1, height: 1 }),
  static convexCapsule_tempVec = vec2.create();
  /**Convex/capsule narrowphase
  * @param {Body}       convexBody
  * @param {Convex}     convexShape
  * @param {Array}      convexPosition
  * @param {Number}     convexAngle
  * @param {Body}       capsuleBody
  * @param {Capsule}    capsuleShape
  * @param {Array}      capsulePosition
  * @param {Number}     capsuleAngle
  */
  convexCapsule(
    convexBody,
    convexShape,
    convexPosition,
    convexAngle,
    capsuleBody,
    capsuleShape,
    capsulePosition,
    capsuleAngle,
    justTest
  ) {
    // Check the circles
    // Add offsets!
    let circlePos = Narrowphase.convexCapsule_tempVec;
    vec2.set(circlePos, capsuleShape.length / 2, 0);
    vec2.rotate(circlePos, circlePos, capsuleAngle);
    vec2.add(circlePos, circlePos, capsulePosition);
    let result1 = this.circleConvex(capsuleBody, capsuleShape, circlePos, capsuleAngle, convexBody, convexShape, convexPosition, convexAngle, justTest, capsuleShape.radius);
    vec2.set(circlePos, -capsuleShape.length / 2, 0);
    vec2.rotate(circlePos, circlePos, capsuleAngle);
    vec2.add(circlePos, circlePos, capsulePosition);
    let result2 = this.circleConvex(capsuleBody, capsuleShape, circlePos, capsuleAngle, convexBody, convexShape, convexPosition, convexAngle, justTest, capsuleShape.radius);
    if (justTest && (result1 || result2)) {
      return true;
    }
    // Check center rect
    let r = Narrowphase.convexCapsule_tempRect;
    Narrowphase.setConvexToCapsuleShapeMiddle(r, capsuleShape);
    let result = this.convexConvex(convexBody, convexShape, convexPosition, convexAngle, capsuleBody, r, capsulePosition, capsuleAngle, justTest);
    return result + result1 + result2;
  }
  /**Capsule/line narrowphase
  * @param {Body}       lineBody
  * @param {Line}       lineShape
  * @param {Array}      linePosition
  * @param {Number}     lineAngle
  * @param {Body}       capsuleBody
  * @param {Capsule}    capsuleShape
  * @param {Array}      capsulePosition
  * @param {Number}     capsuleAngle
  * @todo Implement me!
  */
  lineCapsule(
    lineBody,
    lineShape,
    linePosition,
    lineAngle,
    capsuleBody,
    capsuleShape,
    capsulePosition,
    capsuleAngle,
    justTest
  ) {
    // TODO
    if (justTest) {
      return false;
    } else {
      return 0;
    }
  }
  static capsuleCapsule_tempVec1 = vec2.create();
  static capsuleCapsule_tempVec2 = vec2.create();
  static capsuleCapsule_tempRect1 = new Box({ width: 1, height: 1 });
  /**Capsule/capsule narrowphase
  * @param {Body}       bi
  * @param {Capsule}    si
  * @param {Array}      xi
  * @param {Number}     ai
  * @param {Body}       bj
  * @param {Capsule}    sj
  * @param {Array}      xj
  * @param {Number}     aj
  */
  capsuleCapsule(bi, si, xi, ai, bj, sj, xj, aj, justTest) {
    let enableFrictionBefore;
    // Check the circles
    // Add offsets!
    let circlePosi = Narrowphase.capsuleCapsule_tempVec1,
      circlePosj = capsuleCapsule_tempVec2;
    let numContacts = 0;

    // Need 4 circle checks, between all
    for (let i = 0; i < 2; i++) {
      vec2.set(circlePosi, (i === 0 ? -1 : 1) * si.length / 2, 0);
      vec2.rotate(circlePosi, circlePosi, ai);
      vec2.add(circlePosi, circlePosi, xi);
      for (let j = 0; j < 2; j++) {
        vec2.set(circlePosj, (j === 0 ? -1 : 1) * sj.length / 2, 0);
        vec2.rotate(circlePosj, circlePosj, aj);
        vec2.add(circlePosj, circlePosj, xj);
        // Temporarily turn off friction
        if (this.enableFrictionReduction) {
          enableFrictionBefore = this.enableFriction;
          this.enableFriction = false;
        }
        let result = this.circleCircle(bi, si, circlePosi, ai, bj, sj, circlePosj, aj, justTest, si.radius, sj.radius);
        if (this.enableFrictionReduction) {
          this.enableFriction = enableFrictionBefore;
        }
        if (justTest && result) {
          return true;
        }
        numContacts += result;
      }
    }
    if (this.enableFrictionReduction) {
      // Temporarily turn off friction
      enableFrictionBefore = this.enableFriction;
      this.enableFriction = false;
    }
    // Check circles against the center boxs
    let rect = Narrowphase.capsuleCapsule_tempRect1;
    Narrowphase.setConvexToCapsuleShapeMiddle(rect, si);
    let result1 = this.convexCapsule(bi, rect, xi, ai, bj, sj, xj, aj, justTest);
    if (this.enableFrictionReduction) {
      this.enableFriction = enableFrictionBefore;
    }
    if (justTest && result1) {
      return true;
    }
    numContacts += result1;
    if (this.enableFrictionReduction) {
      // Temporarily turn off friction
      let enableFrictionBefore = this.enableFriction;
      this.enableFriction = false;
    }
    Narrowphase.setConvexToCapsuleShapeMiddle(rect, sj);
    let result2 = this.convexCapsule(bj, rect, xj, aj, bi, si, xi, ai, justTest);
    if (this.enableFrictionReduction) {
      this.enableFriction = enableFrictionBefore;
    }
    if (justTest && result2) {
      return true;
    }
    numContacts += result2;
    if (this.enableFrictionReduction) {
      if (numContacts && this.enableFriction) {
        this.frictionEquations.push(this.createFrictionFromAverage(numContacts));
      }
    }
    return numContacts;
  }
  /**Line/line narrowphase
  * @param {Body}       bodyA
  * @param {Line}       shapeA
  * @param {Array}      positionA
  * @param {Number}     angleA
  * @param {Body}       bodyB
  * @param {Line}       shapeB
  * @param {Array}      positionB
  * @param {Number}     angleB
  * @todo Implement me!
  */
  lineLine(
    bodyA,
    shapeA,
    positionA,
    angleA,
    bodyB,
    shapeB,
    positionB,
    angleB,
    justTest
  ) {
    // TODO
    if (justTest) {
      return false;
    } else {
      return 0;
    }
  }
  /**Plane/line Narrowphase
  * @param {Body}   planeBody
  * @param {Plane}  planeShape
  * @param {Array}  planeOffset
  * @param {Number} planeAngle
  * @param {Body}   lineBody
  * @param {Line}   lineShape
  * @param {Array}  lineOffset
  * @param {Number} lineAngle
  */
  planeLine(planeBody, planeShape, planeOffset, planeAngle,
    lineBody, lineShape, lineOffset, lineAngle, justTest) {
    let worldVertex0 = tmp1,
      worldVertex1 = tmp2,
      worldVertex01 = tmp3,
      worldVertex11 = tmp4,
      worldEdge = tmp5,
      worldEdgeUnit = tmp6,
      dist = tmp7,
      worldNormal = tmp8,
      worldTangent = tmp9,
      verts = tmpArray,
      numContacts = 0;
    // Get start and end points
    vec2.set(worldVertex0, -lineShape.length / 2, 0);
    vec2.set(worldVertex1, lineShape.length / 2, 0);
    // Not sure why we have to use worldVertex*1 here, but it won't work otherwise. Tired.
    vec2.rotate(worldVertex01, worldVertex0, lineAngle);
    vec2.rotate(worldVertex11, worldVertex1, lineAngle);
    add(worldVertex01, worldVertex01, lineOffset);
    add(worldVertex11, worldVertex11, lineOffset);
    vec2.copy(worldVertex0, worldVertex01);
    vec2.copy(worldVertex1, worldVertex11);
    // Get vector along the line
    sub(worldEdge, worldVertex1, worldVertex0);
    vec2.normalize(worldEdgeUnit, worldEdge);
    // Get tangent to the edge.
    vec2.rotate90cw(worldTangent, worldEdgeUnit);
    vec2.rotate(worldNormal, yAxis, planeAngle);
    // Check line ends
    verts[0] = worldVertex0;
    verts[1] = worldVertex1;
    for (let i = 0; i < verts.length; i++) {
      let v = verts[i];
      sub(dist, v, planeOffset);
      let d = dot(dist, worldNormal);
      if (d < 0) {
        if (justTest) {
          return true;
        }
        let c = this.createContactEquation(planeBody, lineBody, planeShape, lineShape);
        numContacts++;
        vec2.copy(c.normalA, worldNormal);
        vec2.normalize(c.normalA, c.normalA);
        // distance vector along plane normal
        vec2.scale(dist, worldNormal, d);
        // Vector from plane center to contact
        sub(c.contactPointA, v, dist);
        sub(c.contactPointA, c.contactPointA, planeBody.position);
        // From line center to contact
        sub(c.contactPointB, v, lineOffset);
        add(c.contactPointB, c.contactPointB, lineOffset);
        sub(c.contactPointB, c.contactPointB, lineBody.position);
        this.contactEquations.push(c);
        if (!this.enableFrictionReduction) {
          if (this.enableFriction) {
            this.frictionEquations.push(this.createFrictionFromContact(c));
          }
        }
      }
    }
    if (justTest) {
      return false;
    }
    if (!this.enableFrictionReduction) {
      if (numContacts && this.enableFriction) {
        this.frictionEquations.push(this.createFrictionFromAverage(numContacts));
      }
    }
    return numContacts;
  }
  particleCapsule(
    particleBody,
    particleShape,
    particlePosition,
    particleAngle,
    capsuleBody,
    capsuleShape,
    capsulePosition,
    capsuleAngle,
    justTest
  ) {
    return this.circleLine(particleBody, particleShape, particlePosition, particleAngle, capsuleBody, capsuleShape, capsulePosition, capsuleAngle, justTest, capsuleShape.radius, 0);
  }
  /**Circle/line Narrowphase
  * @param {Body} circleBody
  * @param {Circle} circleShape
  * @param {Array} circleOffset
  * @param {Number} circleAngle
  * @param {Body} lineBody
  * @param {Line} lineShape
  * @param {Array} lineOffset
  * @param {Number} lineAngle
  * @param {Boolean} justTest If set to true, this function will return the result (intersection or not) without adding equations.
  * @param {Number} lineRadius Radius to add to the line. Can be used to test Capsules.
  * @param {Number} circleRadius If set, this value overrides the circle shape radius.
  */
  circleLine(
    circleBody,
    circleShape,
    circleOffset,
    circleAngle,
    lineBody,
    lineShape,
    lineOffset,
    lineAngle,
    justTest,
    lineRadius,
    circleRadius
  ) {
    let lineRadius = lineRadius || 0,
      circleRadius = typeof (circleRadius) !== "undefined" ? circleRadius : circleShape.radius,
      orthoDist = tmp1,
      lineToCircleOrthoUnit = tmp2,
      projectedPoint = tmp3,
      centerDist = tmp4,
      worldTangent = tmp5,
      worldEdge = tmp6,
      worldEdgeUnit = tmp7,
      worldVertex0 = tmp8,
      worldVertex1 = tmp9,
      worldVertex01 = tmp10,
      worldVertex11 = tmp11,
      dist = tmp12,
      lineToCircle = tmp13,
      lineEndToLineRadius = tmp14,
      verts = tmpArray;
    // Get start and end points
    vec2.set(worldVertex0, -lineShape.length / 2, 0);
    vec2.set(worldVertex1, lineShape.length / 2, 0);
    // Not sure why we have to use worldVertex*1 here, but it won't work otherwise. Tired.
    vec2.rotate(worldVertex01, worldVertex0, lineAngle);
    vec2.rotate(worldVertex11, worldVertex1, lineAngle);
    add(worldVertex01, worldVertex01, lineOffset);
    add(worldVertex11, worldVertex11, lineOffset);
    vec2.copy(worldVertex0, worldVertex01);
    vec2.copy(worldVertex1, worldVertex11);
    // Get vector along the line
    sub(worldEdge, worldVertex1, worldVertex0);
    vec2.normalize(worldEdgeUnit, worldEdge);
    // Get tangent to the edge.
    vec2.rotate90cw(worldTangent, worldEdgeUnit);
    // Check distance from the plane spanned by the edge vs the circle
    sub(dist, circleOffset, worldVertex0);
    let d = dot(dist, worldTangent); // Distance from center of line to circle center
    sub(centerDist, worldVertex0, lineOffset);
    sub(lineToCircle, circleOffset, lineOffset);
    let radiusSum = circleRadius + lineRadius;
    if (Math.abs(d) < radiusSum) {
      // Now project the circle onto the edge
      vec2.scale(orthoDist, worldTangent, d);
      sub(projectedPoint, circleOffset, orthoDist);
      // Add the missing line radius
      vec2.scale(lineToCircleOrthoUnit, worldTangent, dot(worldTangent, lineToCircle));
      vec2.normalize(lineToCircleOrthoUnit, lineToCircleOrthoUnit);
      vec2.scale(lineToCircleOrthoUnit, lineToCircleOrthoUnit, lineRadius);
      add(projectedPoint, projectedPoint, lineToCircleOrthoUnit);
      // Check if the point is within the edge span
      let pos = dot(worldEdgeUnit, projectedPoint);
      let pos0 = dot(worldEdgeUnit, worldVertex0);
      let pos1 = dot(worldEdgeUnit, worldVertex1);
      if (pos > pos0 && pos < pos1) {
        // We got contact!
        if (justTest) {
          return true;
        }
        let c = this.createContactEquation(circleBody, lineBody, circleShape, lineShape);
        vec2.scale(c.normalA, orthoDist, -1);
        vec2.normalize(c.normalA, c.normalA);
        vec2.scale(c.contactPointA, c.normalA, circleRadius);
        add(c.contactPointA, c.contactPointA, circleOffset);
        sub(c.contactPointA, c.contactPointA, circleBody.position);
        sub(c.contactPointB, projectedPoint, lineOffset);
        add(c.contactPointB, c.contactPointB, lineOffset);
        sub(c.contactPointB, c.contactPointB, lineBody.position);
        this.contactEquations.push(c);
        if (this.enableFriction) {
          this.frictionEquations.push(this.createFrictionFromContact(c));
        }
        return 1;
      }
    }
    // Add corner
    verts[0] = worldVertex0;
    verts[1] = worldVertex1;
    for (let i = 0; i < verts.length; i++) {
      let v = verts[i];
      sub(dist, v, circleOffset);
      if (vec2.squaredLength(dist) < Math.pow(radiusSum, 2)) {
        if (justTest) {
          return true;
        }
        let c = this.createContactEquation(circleBody, lineBody, circleShape, lineShape);
        vec2.copy(c.normalA, dist);
        vec2.normalize(c.normalA, c.normalA);
        // Vector from circle to contact point is the normal times the circle radius
        vec2.scale(c.contactPointA, c.normalA, circleRadius);
        add(c.contactPointA, c.contactPointA, circleOffset);
        sub(c.contactPointA, c.contactPointA, circleBody.position);
        sub(c.contactPointB, v, lineOffset);
        vec2.scale(lineEndToLineRadius, c.normalA, -lineRadius);
        add(c.contactPointB, c.contactPointB, lineEndToLineRadius);
        add(c.contactPointB, c.contactPointB, lineOffset);
        sub(c.contactPointB, c.contactPointB, lineBody.position);
        this.contactEquations.push(c);
        if (this.enableFriction) {
          this.frictionEquations.push(this.createFrictionFromContact(c));
        }
        return 1;
      }
    }
    return 0;
  }
  /**Circle/capsule Narrowphase
  * @param {Body}   bi
  * @param {Circle} si
  * @param {Array}  xi
  * @param {Number} ai
  * @param {Body}   bj
  * @param {Line}   sj
  * @param {Array}  xj
  * @param {Number} aj
  */
  circleCapsule(bi, si, xi, ai, bj, sj, xj, aj, justTest) {
    return this.circleLine(bi, si, xi, ai, bj, sj, xj, aj, justTest, sj.radius);
  }
  /**Circle/convex Narrowphase.
  * @param {Body} circleBody
  * @param {Circle} circleShape
  * @param {Array} circleOffset
  * @param {Number} circleAngle
  * @param {Body} convexBody
  * @param {Convex} convexShape
  * @param {Array} convexOffset
  * @param {Number} convexAngle
  * @param {Boolean} justTest
  * @param {Number} circleRadius
  */
  circleConvex(
    circleBody,
    circleShape,
    circleOffset,
    circleAngle,
    convexBody,
    convexShape,
    convexOffset,
    convexAngle,
    justTest,
    circleRadius
  ) {
    let circleRadius = typeof (circleRadius) === "number" ? circleRadius : circleShape.radius;
    let worldVertex0 = tmp1,
      worldVertex1 = tmp2,
      worldEdge = tmp3,
      worldEdgeUnit = tmp4,
      worldNormal = tmp5,
      centerDist = tmp6,
      convexToCircle = tmp7,
      orthoDist = tmp8,
      projectedPoint = tmp9,
      dist = tmp10,
      worldVertex = tmp11,
      closestEdge = -1,
      closestEdgeDistance = null,
      closestEdgeOrthoDist = tmp12,
      closestEdgeProjectedPoint = tmp13,
      candidate = tmp14,
      candidateDist = tmp15,
      minCandidate = tmp16,
      found = false,
      minCandidateDistance = Number.MAX_VALUE;
    let numReported = 0;
    // New algorithm:
    // 1. Check so center of circle is not inside the polygon. If it is, this wont work...
    // 2. For each edge
    // 2. 1. Get point on circle that is closest to the edge (scale normal with -radius)
    // 2. 2. Check if point is inside.
    let verts = convexShape.vertices;
    // Check all edges first
    for (let i = 0; i !== verts.length + 1; i++) {
      let Ray.v0 = verts[i % verts.length],
        v1 = verts[(i + 1) % verts.length];
      vec2.rotate(worldVertex0, Ray.v0, convexAngle);
      vec2.rotate(worldVertex1, v1, convexAngle);
      add(worldVertex0, worldVertex0, convexOffset);
      add(worldVertex1, worldVertex1, convexOffset);
      sub(worldEdge, worldVertex1, worldVertex0);
      vec2.normalize(worldEdgeUnit, worldEdge);
      // Get tangent to the edge. Points out of the Convex
      vec2.rotate90cw(worldNormal, worldEdgeUnit);
      // Get point on circle, closest to the polygon
      vec2.scale(candidate, worldNormal, -circleShape.radius);
      add(candidate, candidate, circleOffset);
      if (pointInConvex(candidate, convexShape, convexOffset, convexAngle)) {
        vec2.sub(candidateDist, worldVertex0, candidate);
        let candidateDistance = Math.abs(vec2.dot(candidateDist, worldNormal));
        if (candidateDistance < minCandidateDistance) {
          vec2.copy(minCandidate, candidate);
          minCandidateDistance = candidateDistance;
          vec2.scale(closestEdgeProjectedPoint, worldNormal, candidateDistance);
          vec2.add(closestEdgeProjectedPoint, closestEdgeProjectedPoint, candidate);
          found = true;
        }
      }
    }
    if (found) {
      if (justTest) {
        return true;
      }
      let c = this.createContactEquation(circleBody, convexBody, circleShape, convexShape);
      vec2.sub(c.normalA, minCandidate, circleOffset);
      vec2.normalize(c.normalA, c.normalA);
      vec2.scale(c.contactPointA, c.normalA, circleRadius);
      add(c.contactPointA, c.contactPointA, circleOffset);
      sub(c.contactPointA, c.contactPointA, circleBody.position);
      sub(c.contactPointB, closestEdgeProjectedPoint, convexOffset);
      add(c.contactPointB, c.contactPointB, convexOffset);
      sub(c.contactPointB, c.contactPointB, convexBody.position);
      this.contactEquations.push(c);
      if (this.enableFriction) {
        this.frictionEquations.push(this.createFrictionFromContact(c));
      }
      return 1;
    }
    // Check all vertices
    if (circleRadius > 0) {
      for (let i = 0; i < verts.length; i++) {
        let localVertex = verts[i];
        vec2.rotate(worldVertex, localVertex, convexAngle);
        add(worldVertex, worldVertex, convexOffset);
        sub(dist, worldVertex, circleOffset);
        if (vec2.squaredLength(dist) < Math.pow(circleRadius, 2)) {
          if (justTest) {
            return true;
          }
          let c = this.createContactEquation(circleBody, convexBody, circleShape, convexShape);
          vec2.copy(c.normalA, dist);
          vec2.normalize(c.normalA, c.normalA);
          // Vector from circle to contact point is the normal times the circle radius
          vec2.scale(c.contactPointA, c.normalA, circleRadius);
          add(c.contactPointA, c.contactPointA, circleOffset);
          sub(c.contactPointA, c.contactPointA, circleBody.position);
          sub(c.contactPointB, worldVertex, convexOffset);
          add(c.contactPointB, c.contactPointB, convexOffset);
          sub(c.contactPointB, c.contactPointB, convexBody.position);
          this.contactEquations.push(c);
          if (this.enableFriction) {
            this.frictionEquations.push(this.createFrictionFromContact(c));
          }
          return 1;
        }
      }
    }
    return 0;
  }
  static pic_worldVertex0 = vec2.create();
  static pic_worldVertex1 = vec2.create();
  static pic_r0 = vec2.create();
  static pic_r1 = vec2.create();
  /*Check if a point is in a polygon
  */
  static pointInConvex(worldPoint, convexShape, convexOffset, convexAngle) {
    let worldVertex0 = Narrowphase.pic_worldVertex0,
      worldVertex1 = Narrowphase.pic_worldVertex1,
      r0 = Narrowphase.pic_r0,
      r1 = Narrowphase.pic_r1,
      point = worldPoint,
      verts = convexShape.vertices,
      lastCross = null;
    for (let i = 0; i !== verts.length + 1; i++) {
      let Ray.v0 = verts[i % verts.length],
        v1 = verts[(i + 1) % verts.length];
      // Transform vertices to world
      // @todo The point should be transformed to local coordinates in the convex, no need to transform each vertex
      vec2.rotate(worldVertex0, Ray.v0, convexAngle);
      vec2.rotate(worldVertex1, v1, convexAngle);
      add(worldVertex0, worldVertex0, convexOffset);
      add(worldVertex1, worldVertex1, convexOffset);
      sub(r0, worldVertex0, point);
      sub(r1, worldVertex1, point);
      let cross = vec2.crossLength(r0, r1);
      if (lastCross === null) {
        lastCross = cross;
      }
      // If we got a different sign of the distance vector, the point is out of the polygon
      if (cross * lastCross <= 0) {
        return false;
      }
      lastCross = cross;
    }
    return true;
  }
  /**Particle/convex Narrowphase
  * @param {Body} particleBody
  * @param {Particle} particleShape
  * @param {Array} particleOffset
  * @param {Number} particleAngle
  * @param {Body} convexBody
  * @param {Convex} convexShape
  * @param {Array} convexOffset
  * @param {Number} convexAngle
  * @param {Boolean} justTest
  * @todo use pointInConvex and code more similar to circleConvex
  * @todo don't transform each vertex, but transform the particle position to convex-local instead
  */
  particleConvex(
    particleBody,
    particleShape,
    particleOffset,
    particleAngle,
    convexBody,
    convexShape,
    convexOffset,
    convexAngle,
    justTest
  ) {
    let worldVertex0 = tmp1,
      worldVertex1 = tmp2,
      worldEdge = tmp3,
      worldEdgeUnit = tmp4,
      worldTangent = tmp5,
      centerDist = tmp6,
      convexToparticle = tmp7,
      orthoDist = tmp8,
      projectedPoint = tmp9,
      dist = tmp10,
      worldVertex = tmp11,
      closestEdge = -1,
      closestEdgeDistance = null,
      closestEdgeOrthoDist = tmp12,
      closestEdgeProjectedPoint = tmp13,
      r0 = tmp14, // vector from particle to vertex0
      r1 = tmp15,
      localPoint = tmp16,
      candidateDist = tmp17,
      minEdgeNormal = tmp18,
      minCandidateDistance = Number.MAX_VALUE;
    let numReported = 0,
      found = false,
      verts = convexShape.vertices;
    // Check if the particle is in the polygon at all
    if (!pointInConvex(particleOffset, convexShape, convexOffset, convexAngle)) {
      return 0;
    }
    if (justTest) {
      return true;
    }
    // Check edges first
    let lastCross = null;
    for (let i = 0; i !== verts.length + 1; i++) {
      let Ray.v0 = verts[i % verts.length],
        v1 = verts[(i + 1) % verts.length];
      // Transform vertices to world
      vec2.rotate(worldVertex0, Ray.v0, convexAngle);
      vec2.rotate(worldVertex1, v1, convexAngle);
      add(worldVertex0, worldVertex0, convexOffset);
      add(worldVertex1, worldVertex1, convexOffset);
      // Get world edge
      sub(worldEdge, worldVertex1, worldVertex0);
      vec2.normalize(worldEdgeUnit, worldEdge);
      // Get tangent to the edge. Points out of the Convex
      vec2.rotate90cw(worldTangent, worldEdgeUnit);
      // Check distance from the infinite line (spanned by the edge) to the particle
      sub(dist, particleOffset, worldVertex0);
      let d = dot(dist, worldTangent);
      sub(centerDist, worldVertex0, convexOffset);
      sub(convexToparticle, particleOffset, convexOffset);
      vec2.sub(candidateDist, worldVertex0, particleOffset);
      let candidateDistance = Math.abs(vec2.dot(candidateDist, worldTangent));
      if (candidateDistance < minCandidateDistance) {
        minCandidateDistance = candidateDistance;
        vec2.scale(closestEdgeProjectedPoint, worldTangent, candidateDistance);
        vec2.add(closestEdgeProjectedPoint, closestEdgeProjectedPoint, particleOffset);
        vec2.copy(minEdgeNormal, worldTangent);
        found = true;
      }
    }
    if (found) {
      let c = this.createContactEquation(particleBody, convexBody, particleShape, convexShape);
      vec2.scale(c.normalA, minEdgeNormal, -1);
      vec2.normalize(c.normalA, c.normalA);
      // Particle has no extent to the contact point
      vec2.set(c.contactPointA, 0, 0);
      add(c.contactPointA, c.contactPointA, particleOffset);
      sub(c.contactPointA, c.contactPointA, particleBody.position);
      // From convex center to point
      sub(c.contactPointB, closestEdgeProjectedPoint, convexOffset);
      add(c.contactPointB, c.contactPointB, convexOffset);
      sub(c.contactPointB, c.contactPointB, convexBody.position);
      this.contactEquations.push(c);
      if (this.enableFriction) {
        this.frictionEquations.push(this.createFrictionFromContact(c));
      }
      return 1;
    }

    return 0;
  }
  /**Circle/circle Narrowphase
  * @param {Body} bodyA
  * @param {Circle} shapeA
  * @param {Array} offsetA
  * @param {Number} angleA
  * @param {Body} bodyB
  * @param {Circle} shapeB
  * @param {Array} offsetB
  * @param {Number} angleB
  * @param {Boolean} justTest
  * @param {Number} [radiusA] Optional radius to use for shapeA
  * @param {Number} [radiusB] Optional radius to use for shapeB
  */
  circleCircle(
    bodyA,
    shapeA,
    offsetA,
    angleA,
    bodyB,
    shapeB,
    offsetB,
    angleB,
    justTest,
    radiusA,
    radiusB
  ) {
    let dist = tmp1,
      radiusA = radiusA || shapeA.radius,
      radiusB = radiusB || shapeB.radius;
    sub(dist, offsetA, offsetB);
    let r = radiusA + radiusB;
    if (vec2.squaredLength(dist) > Math.pow(r, 2)) {
      return 0;
    }
    if (justTest) {
      return true;
    }
    let c = this.createContactEquation(bodyA, bodyB, shapeA, shapeB);
    sub(c.normalA, offsetB, offsetA);
    vec2.normalize(c.normalA, c.normalA);
    vec2.scale(c.contactPointA, c.normalA, radiusA);
    vec2.scale(c.contactPointB, c.normalA, -radiusB);
    add(c.contactPointA, c.contactPointA, offsetA);
    sub(c.contactPointA, c.contactPointA, bodyA.position);
    add(c.contactPointB, c.contactPointB, offsetB);
    sub(c.contactPointB, c.contactPointB, bodyB.position);
    this.contactEquations.push(c);
    if (this.enableFriction) {
      this.frictionEquations.push(this.createFrictionFromContact(c));
    }
    return 1;
  }
  /**Plane/Convex Narrowphase
  * @param {Body} planeBody
  * @param {Plane} planeShape
  * @param {Array} planeOffset
  * @param {Number} planeAngle
  * @param {Body} convexBody
  * @param {Convex} convexShape
  * @param {Array} convexOffset
  * @param {Number} convexAngle
  * @param {Boolean} justTest
  */
  planeConvex(
    planeBody,
    planeShape,
    planeOffset,
    planeAngle,
    convexBody,
    convexShape,
    convexOffset,
    convexAngle,
    justTest
  ) {
    let worldVertex = tmp1,
      worldNormal = tmp2,
      dist = tmp3;
    let numReported = 0;
    vec2.rotate(worldNormal, yAxis, planeAngle);
    for (let i = 0; i !== convexShape.vertices.length; i++) {
      let v = convexShape.vertices[i];
      vec2.rotate(worldVertex, v, convexAngle);
      add(worldVertex, worldVertex, convexOffset);
      sub(dist, worldVertex, planeOffset);
      if (dot(dist, worldNormal) <= 0) {
        if (justTest) {
          return true;
        }
        // Found vertex
        numReported++;
        let c = this.createContactEquation(planeBody, convexBody, planeShape, convexShape);
        sub(dist, worldVertex, planeOffset);
        vec2.copy(c.normalA, worldNormal);
        let d = dot(dist, c.normalA);
        vec2.scale(dist, c.normalA, d);
        // rj is from convex center to contact
        sub(c.contactPointB, worldVertex, convexBody.position);

        // ri is from plane center to contact
        sub(c.contactPointA, worldVertex, dist);
        sub(c.contactPointA, c.contactPointA, planeBody.position);
        this.contactEquations.push(c);
        if (!this.enableFrictionReduction) {
          if (this.enableFriction) {
            this.frictionEquations.push(this.createFrictionFromContact(c));
          }
        }
      }
    }
    if (this.enableFrictionReduction) {
      if (this.enableFriction && numReported) {
        this.frictionEquations.push(this.createFrictionFromAverage(numReported));
      }
    }
    return numReported;
  }
  /**Narrowphase for particle vs plane
  * @param {Body}       particleBody
  * @param {Particle}   particleShape
  * @param {Array}      particleOffset
  * @param {Number}     particleAngle
  * @param {Body}       planeBody
  * @param {Plane}      planeShape
  * @param {Array}      planeOffset
  * @param {Number}     planeAngle
  * @param {Boolean}     justTest
  */
  particlePlane(
    particleBody,
    particleShape,
    particleOffset,
    particleAngle,
    planeBody,
    planeShape,
    planeOffset,
    planeAngle,
    justTest
  ) {
    let dist = tmp1,
      worldNormal = tmp2;
    planeAngle = planeAngle || 0;
    sub(dist, particleOffset, planeOffset);
    vec2.rotate(worldNormal, yAxis, planeAngle);
    let d = dot(dist, worldNormal);
    if (d > 0) {
      return 0;
    }
    if (justTest) {
      return true;
    }
    let c = this.createContactEquation(planeBody, particleBody, planeShape, particleShape);
    vec2.copy(c.normalA, worldNormal);
    vec2.scale(dist, c.normalA, d);
    // dist is now the distance vector in the normal direction
    // ri is the particle position projected down onto the plane, from the plane center
    sub(c.contactPointA, particleOffset, dist);
    sub(c.contactPointA, c.contactPointA, planeBody.position);
    // rj is from the body center to the particle center
    sub(c.contactPointB, particleOffset, particleBody.position);
    this.contactEquations.push(c);
    if (this.enableFriction) {
      this.frictionEquations.push(this.createFrictionFromContact(c));
    }
    return 1;
  }
  /**Circle/Particle Narrowphase
  * @param {Body} circleBody
  * @param {Circle} circleShape
  * @param {Array} circleOffset
  * @param {Number} circleAngle
  * @param {Body} particleBody
  * @param {Particle} particleShape
  * @param {Array} particleOffset
  * @param {Number} particleAngle
  * @param {Boolean} justTest
  */
  circleParticle(
    circleBody,
    circleShape,
    circleOffset,
    circleAngle,
    particleBody,
    particleShape,
    particleOffset,
    particleAngle,
    justTest
  ) {
    let dist = tmp1;
    sub(dist, particleOffset, circleOffset);
    if (vec2.squaredLength(dist) > Math.pow(circleShape.radius, 2)) {
      return 0;
    }
    if (justTest) {
      return true;
    }
    let c = this.createContactEquation(circleBody, particleBody, circleShape, particleShape);
    vec2.copy(c.normalA, dist);
    vec2.normalize(c.normalA, c.normalA);
    // Vector from circle to contact point is the normal times the circle radius
    vec2.scale(c.contactPointA, c.normalA, circleShape.radius);
    add(c.contactPointA, c.contactPointA, circleOffset);
    sub(c.contactPointA, c.contactPointA, circleBody.position);
    // Vector from particle center to contact point is zero
    sub(c.contactPointB, particleOffset, particleBody.position);
    this.contactEquations.push(c);
    if (this.enableFriction) {
      this.frictionEquations.push(this.createFrictionFromContact(c));
    }
    return 1;
  }
  static planeCapsule_tmpCircle = new Circle({ radius: 1 });
  static planeCapsule_tmp1 = vec2.create();
  static planeCapsule_tmp2 = vec2.create();
  static planeCapsule_tmp3 = vec2.create();
  /**@param {Body} planeBody
  * @param {Circle} planeShape
  * @param {Array} planeOffset
  * @param {Number} planeAngle
  * @param {Body} capsuleBody
  * @param {Particle} capsuleShape
  * @param {Array} capsuleOffset
  * @param {Number} capsuleAngle
  * @param {Boolean} justTest
  */
  planeCapsule(
    planeBody,
    planeShape,
    planeOffset,
    planeAngle,
    capsuleBody,
    capsuleShape,
    capsuleOffset,
    capsuleAngle,
    justTest
  ) {
    let end1 = Narrowphase.planeCapsule_tmp1,
      end2 = Narrowphase.planeCapsule_tmp2,
      circle = Narrowphase.planeCapsule_tmpCircle,
      dst = Narrowphase.planeCapsule_tmp3;
    // Compute world end positions
    vec2.set(end1, -capsuleShape.length / 2, 0);
    vec2.rotate(end1, end1, capsuleAngle);
    add(end1, end1, capsuleOffset);
    vec2.set(end2, capsuleShape.length / 2, 0);
    vec2.rotate(end2, end2, capsuleAngle);
    add(end2, end2, capsuleOffset);
    circle.radius = capsuleShape.radius;
    let enableFrictionBefore;
    // Temporarily turn off friction
    if (this.enableFrictionReduction) {
      enableFrictionBefore = this.enableFriction;
      this.enableFriction = false;
    }
    // Do Narrowphase as two circles
    let numContacts1 = this.circlePlane(capsuleBody, circle, end1, 0, planeBody, planeShape, planeOffset, planeAngle, justTest),
      numContacts2 = this.circlePlane(capsuleBody, circle, end2, 0, planeBody, planeShape, planeOffset, planeAngle, justTest);
    // Restore friction
    if (this.enableFrictionReduction) {
      this.enableFriction = enableFrictionBefore;
    }
    if (justTest) {
      return numContacts1 || numContacts2;
    } else {
      let numTotal = numContacts1 + numContacts2;
      if (this.enableFrictionReduction) {
        if (numTotal) {
          this.frictionEquations.push(this.createFrictionFromAverage(numTotal));
        }
      }
      return numTotal;
    }
  }
  /**Creates ContactEquations and FrictionEquations for a collision.
  * @param {Body}    bi     The first body that should be connected to the equations.
  * @param {Circle}  si     The circle shape participating in the collision.
  * @param {Array}   xi     Extra offset to take into account for the Shape, in addition to the one in circleBody.position. Will *not* be rotated by circleBody.angle (maybe it should, for sake of homogenity?). Set to null if none.
  * @param {Body}    bj     The second body that should be connected to the equations.
  * @param {Plane}   sj     The Plane shape that is participating
  * @param {Array}   xj     Extra offset for the plane shape.
  * @param {Number}  aj     Extra angle to apply to the plane
  */
  circlePlane(bi, si, xi, ai, bj, sj, xj, aj, justTest) {
    let circleBody = bi,
      circleShape = si,
      circleOffset = xi, // Offset from body center, rotated!
      planeBody = bj,
      shapeB = sj,
      planeOffset = xj,
      planeAngle = aj;
    planeAngle = planeAngle || 0;
    // Vector from plane to circle
    let planeToCircle = tmp1,
      worldNormal = tmp2,
      temp = tmp3;
    sub(planeToCircle, circleOffset, planeOffset);
    // World plane normal
    vec2.rotate(worldNormal, yAxis, planeAngle);
    // Normal direction distance
    let d = dot(worldNormal, planeToCircle);
    if (d > circleShape.radius) {
      return 0; // No overlap. Abort.
    }
    if (justTest) {
      return true;
    }
    // Create contact
    let contact = this.createContactEquation(planeBody, circleBody, sj, si);
    // ni is the plane world normal
    vec2.copy(contact.normalA, worldNormal);
    // rj is the vector from circle center to the contact point
    vec2.scale(contact.contactPointB, contact.normalA, -circleShape.radius);
    add(contact.contactPointB, contact.contactPointB, circleOffset);
    sub(contact.contactPointB, contact.contactPointB, circleBody.position);
    // ri is the distance from plane center to contact.
    vec2.scale(temp, contact.normalA, d);
    sub(contact.contactPointA, planeToCircle, temp); // Subtract normal distance vector from the distance vector
    add(contact.contactPointA, contact.contactPointA, planeOffset);
    sub(contact.contactPointA, contact.contactPointA, planeBody.position);
    this.contactEquations.push(contact);
    if (this.enableFriction) {
      this.frictionEquations.push(this.createFrictionFromContact(contact));
    }
    return 1;
  }
  /**Convex/convex Narrowphase.See <a href="http://www.altdevblogaday.com/2011/05/13/contact-generation-between-3d-convex-meshes/">this article</a> for more info.
  * @param {Body} bi
  * @param {Convex} si
  * @param {Array} xi
  * @param {Number} ai
  * @param {Body} bj
  * @param {Convex} sj
  * @param {Array} xj
  * @param {Number} aj
  */
  convexConvex(bi, si, xi, ai, bj, sj, xj, aj, justTest, precision) {
    let sepAxis = tmp1,
      worldPoint = tmp2,
      worldPoint0 = tmp3,
      worldPoint1 = tmp4,
      worldEdge = tmp5,
      projected = tmp6,
      penetrationVec = tmp7,
      dist = tmp8,
      worldNormal = tmp9,
      numContacts = 0,
      precision = typeof (precision) === 'number' ? precision : 0;
    let found = Narrowphase.findSeparatingAxis(si, xi, ai, sj, xj, aj, sepAxis);
    if (!found) {
      return 0;
    }
    // Make sure the separating axis is directed from shape i to shape j
    sub(dist, xj, xi);
    if (dot(sepAxis, dist) > 0) {
      vec2.scale(sepAxis, sepAxis, -1);
    }
    // Find edges with normals closest to the separating axis
    let closestEdge1 = Narrowphase.getClosestEdge(si, ai, sepAxis, true), // Flipped axis
      closestEdge2 = Narrowphase.getClosestEdge(sj, aj, sepAxis);
    if (closestEdge1 === -1 || closestEdge2 === -1) {
      return 0;
    }
    // Loop over the shapes
    for (let k = 0; k < 2; k++) {
      let closestEdgeA = closestEdge1,
        closestEdgeB = closestEdge2,
        shapeA = si, shapeB = sj,
        offsetA = xi, offsetB = xj,
        angleA = ai, angleB = aj,
        bodyA = bi, bodyB = bj;
      if (k === 0) {
        // Swap!
        let tmp;
        tmp = closestEdgeA;
        closestEdgeA = closestEdgeB;
        closestEdgeB = tmp;
        tmp = shapeA;
        shapeA = shapeB;
        shapeB = tmp;
        tmp = offsetA;
        offsetA = offsetB;
        offsetB = tmp;
        tmp = angleA;
        angleA = angleB;
        angleB = tmp;
        tmp = bodyA;
        bodyA = bodyB;
        bodyB = tmp;
      }
      // Loop over 2 points in convex B
      for (let j = closestEdgeB; j < closestEdgeB + 2; j++) {
        // Get world point
        let v = shapeB.vertices[(j + shapeB.vertices.length) % shapeB.vertices.length];
        vec2.rotate(worldPoint, v, angleB);
        add(worldPoint, worldPoint, offsetB);
        let insideNumEdges = 0;
        // Loop over the 3 closest edges in convex A
        for (let i = closestEdgeA - 1; i < closestEdgeA + 2; i++) {
          let Ray.v0 = shapeA.vertices[(i + shapeA.vertices.length) % shapeA.vertices.length],
            v1 = shapeA.vertices[(i + 1 + shapeA.vertices.length) % shapeA.vertices.length];
          // Construct the edge
          vec2.rotate(worldPoint0, Ray.v0, angleA);
          vec2.rotate(worldPoint1, v1, angleA);
          add(worldPoint0, worldPoint0, offsetA);
          add(worldPoint1, worldPoint1, offsetA);
          sub(worldEdge, worldPoint1, worldPoint0);
          vec2.rotate90cw(worldNormal, worldEdge); // Normal points out of convex 1
          vec2.normalize(worldNormal, worldNormal);
          sub(dist, worldPoint, worldPoint0);
          let d = dot(worldNormal, dist);
          if ((i === closestEdgeA && d <= precision) || (i !== closestEdgeA && d <= 0)) {
            insideNumEdges++;
          }
        }
        if (insideNumEdges >= 3) {
          if (justTest) {
            return true;
          }
          // worldPoint was on the "inside" side of each of the 3 checked edges.
          // Project it to the center edge and use the projection direction as normal
          // Create contact
          let c = this.createContactEquation(bodyA, bodyB, shapeA, shapeB);
          numContacts++;
          // Get center edge from body A
          let Ray.v0 = shapeA.vertices[(closestEdgeA) % shapeA.vertices.length],
            v1 = shapeA.vertices[(closestEdgeA + 1) % shapeA.vertices.length];
          // Construct the edge
          vec2.rotate(worldPoint0, Ray.v0, angleA);
          vec2.rotate(worldPoint1, v1, angleA);
          add(worldPoint0, worldPoint0, offsetA);
          add(worldPoint1, worldPoint1, offsetA);
          sub(worldEdge, worldPoint1, worldPoint0);
          vec2.rotate90cw(c.normalA, worldEdge); // Normal points out of convex A
          vec2.normalize(c.normalA, c.normalA);
          sub(dist, worldPoint, worldPoint0); // From edge point to the penetrating point
          let d = dot(c.normalA, dist);             // Penetration
          vec2.scale(penetrationVec, c.normalA, d);     // Vector penetration
          sub(c.contactPointA, worldPoint, offsetA);
          sub(c.contactPointA, c.contactPointA, penetrationVec);
          add(c.contactPointA, c.contactPointA, offsetA);
          sub(c.contactPointA, c.contactPointA, bodyA.position);
          sub(c.contactPointB, worldPoint, offsetB);
          add(c.contactPointB, c.contactPointB, offsetB);
          sub(c.contactPointB, c.contactPointB, bodyB.position);
          this.contactEquations.push(c);
          // Todo reduce to 1 friction equation if we have 2 contact points
          if (!this.enableFrictionReduction) {
            if (this.enableFriction) {
              this.frictionEquations.push(this.createFrictionFromContact(c));
            }
          }
        }
      }
    }
    if (this.enableFrictionReduction) {
      if (this.enableFriction && numContacts) {
        this.frictionEquations.push(this.createFrictionFromAverage(numContacts));
      }
    }
    return numContacts;
  }
  // .projectConvex is called by other functions, need local tmp vectors
  static pcoa_tmp1 = vec2.fromValues(0, 0);
  /**Project a Convex onto a world-oriented axis
    * @param {Convex} convexShape
    * @param {Array} convexOffset
    * @param {Number} convexAngle
    * @param {Array} worldAxis
    * @param {Array} result
    */
  static projectConvexOntoAxis(convexShape, convexOffset, convexAngle, worldAxis, result) {
    let max = null,
      min = null,
      v,
      value,
      localAxis = Narrowphase.pcoa_tmp1;
    // Convert the axis to local coords of the body
    vec2.rotate(localAxis, worldAxis, -convexAngle);
    // Get projected position of all vertices
    for (let i = 0; i < convexShape.vertices.length; i++) {
      v = convexShape.vertices[i];
      value = dot(v, localAxis);
      if (max === null || value > max) {
        max = value;
      }
      if (min === null || value < min) {
        min = value;
      }
    }
    if (min > max) {
      let t = min;
      min = max;
      max = t;
    }
    // Project the position of the body onto the axis - need to add this to the result
    let offset = dot(convexOffset, worldAxis);
    vec2.set(result, min + offset, max + offset);
  }
  // .findSeparatingAxis is called by other functions, need local tmp vectors
  static fsa_tmp1 = vec2.fromValues(0, 0);
  static fsa_tmp2 = vec2.fromValues(0, 0);
  static fsa_tmp3 = vec2.fromValues(0, 0);
  static fsa_tmp4 = vec2.fromValues(0, 0);
  static fsa_tmp5 = vec2.fromValues(0, 0);
  static fsa_tmp6 = vec2.fromValues(0, 0);
  /**Find a separating axis between the shapes, that maximizes the separating distance between them.
  * @param {Convex}     c1
  * @param {Array}      offset1
  * @param {Number}     angle1
  * @param {Convex}     c2
  * @param {Array}      offset2
  * @param {Number}     angle2
  * @param {Array}      sepAxis     The resulting axis
  * @returns {Boolean}                Whether the axis could be found.
  */
  static findSeparatingAxis(c1, offset1, angle1, c2, offset2, angle2, sepAxis) {
    let maxDist = null,
      overlap = false,
      found = false,
      edge = Narrowphase.fsa_tmp1,
      worldPoint0 = Narrowphase.fsa_tmp2,
      worldPoint1 = Narrowphase.fsa_tmp3,
      normal = Narrowphase.fsa_tmp4,
      span1 = Narrowphase.fsa_tmp5,
      span2 = Narrowphase.fsa_tmp6;
    if (c1 instanceof Box && c2 instanceof Box) {
      for (let j = 0; j !== 2; j++) {
        let c = c1,
          angle = angle1;
        if (j === 1) {
          c = c2;
          angle = angle2;
        }
        for (let i = 0; i !== 2; i++) {
          // Get the world edge
          if (i === 0) {
            vec2.set(normal, 0, 1);
          } else if (i === 1) {
            vec2.set(normal, 1, 0);
          }
          if (angle !== 0) {
            vec2.rotate(normal, normal, angle);
          }
          // Project hulls onto that normal
          Narrowphase.projectConvexOntoAxis(c1, offset1, angle1, normal, span1);
          Narrowphase.projectConvexOntoAxis(c2, offset2, angle2, normal, span2);
          // Order by span position
          let a = span1,
            b = span2,
            swapped = false;
          if (span1[0] > span2[0]) {
            b = span1;
            a = span2;
            swapped = true;
          }
          // Get separating distance
          let dist = b[0] - a[1];
          overlap = (dist <= 0);
          if (maxDist === null || dist > maxDist) {
            vec2.copy(sepAxis, normal);
            maxDist = dist;
            found = overlap;
          }
        }
      }
    } else {
      for (let j = 0; j !== 2; j++) {
        let c = c1,
          angle = angle1;
        if (j === 1) {
          c = c2;
          angle = angle2;
        }
        for (let i = 0; i !== c.vertices.length; i++) {
          // Get the world edge
          vec2.rotate(worldPoint0, c.vertices[i], angle);
          vec2.rotate(worldPoint1, c.vertices[(i + 1) % c.vertices.length], angle);
          sub(edge, worldPoint1, worldPoint0);
          // Get normal - just rotate 90 degrees since vertices are given in CCW
          vec2.rotate90cw(normal, edge);
          vec2.normalize(normal, normal);
          // Project hulls onto that normal
          Narrowphase.projectConvexOntoAxis(c1, offset1, angle1, normal, span1);
          Narrowphase.projectConvexOntoAxis(c2, offset2, angle2, normal, span2);
          // Order by span position
          let a = span1,
            b = span2,
            swapped = false;
          if (span1[0] > span2[0]) {
            b = span1;
            a = span2;
            swapped = true;
          }
          // Get separating distance
          let dist = b[0] - a[1];
          overlap = (dist <= 0);
          if (maxDist === null || dist > maxDist) {
            vec2.copy(sepAxis, normal);
            maxDist = dist;
            found = overlap;
          }
        }
      }
    }

    /*
    // Needs to be tested some more
    for(let j=0; j!==2; j++){
        let c = c1,
            angle = angle1;
        if(j===1){
            c = c2;
            angle = angle2;
        }
        for(let i=0; i!==c.axes.length; i++){
            let normal = c.axes[i];
            // Project hulls onto that normal
            Narrowphase.projectConvexOntoAxis(c1, offset1, angle1, normal, span1);
            Narrowphase.projectConvexOntoAxis(c2, offset2, angle2, normal, span2);
            // Order by span position
            let a=span1,
                b=span2,
                swapped = false;
            if(span1[0] > span2[0]){
                b=span1;
                a=span2;
                swapped = true;
            }
            // Get separating distance
            let dist = b[0] - a[1];
            overlap = (dist <= Narrowphase.convexPrecision);
            if(maxDist===null || dist > maxDist){
                vec2.copy(sepAxis, normal);
                maxDist = dist;
                found = overlap;
            }
        }
    }
    */
    return found;
  }
  // .getClosestEdge is called by other functions, need local tmp vectors
  static gce_tmp1 = vec2.fromValues(0, 0);
  static gce_tmp2 = vec2.fromValues(0, 0);
  static gce_tmp3 = vec2.fromValues(0, 0);
  /**Get the edge that has a normal closest to an axis.
  * @param {Convex}     c
  * @param {Number}     angle
  * @param {Array}      axis
  * @param {Boolean}    flip
  * @returns {Number}             Index of the edge that is closest. This index and the next spans the resulting edge. Returns -1 if failed.
  */
  getClosestEdge(c, angle, axis, flip) {
    let localAxis = Narrowphase.gce_tmp1,
      edge = Narrowphase.gce_tmp2,
      normal = Narrowphase.gce_tmp3;
    // Convert the axis to local coords of the body
    vec2.rotate(localAxis, axis, -angle);
    if (flip) {
      vec2.scale(localAxis, localAxis, -1);
    }
    let closestEdge = -1,
      N = c.vertices.length,
      maxDot = -1;
    for (let i = 0; i !== N; i++) {
      // Get the edge
      sub(edge, c.vertices[(i + 1) % N], c.vertices[i % N]);
      // Get normal - just rotate 90 degrees since vertices are given in CCW
      vec2.rotate90cw(normal, edge);
      vec2.normalize(normal, normal);
      let d = dot(normal, localAxis);
      if (closestEdge === -1 || d > maxDot) {
        closestEdge = i % N;
        maxDot = d;
      }
    }
    return closestEdge;
  }
  static circleHeightfield_candidate = vec2.create();
  static circleHeightfield_dist = vec2.create();
  static circleHeightfield_Ray.v0 = vec2.create();
  static circleHeightfield_v1 = vec2.create();
  static circleHeightfield_minCandidate = vec2.create();
  static circleHeightfield_worldNormal = vec2.create();
  static circleHeightfield_minCandidateNormal = vec2.create();
  /**@param {Body}           bi
  * @param {Circle}         si
  * @param {Array}          xi
  * @param {Body}           bj
  * @param {Heightfield}    sj
  * @param {Array}          xj
  * @param {Number}         aj
  */
  circleHeightfield(circleBody, circleShape, circlePos, circleAngle,
    hfBody, hfShape, hfPos, hfAngle, justTest, radius) {
    let data = hfShape.heights,
      radius = radius || circleShape.radius,
      w = hfShape.elementWidth,
      dist = Narrowphase.circleHeightfield_dist,
      candidate = Narrowphase.circleHeightfield_candidate,
      minCandidate = Narrowphase.circleHeightfield_minCandidate,
      minCandidateNormal = Narrowphase.circleHeightfield_minCandidateNormal,
      worldNormal = Narrowphase.circleHeightfield_worldNormal,
      Ray.v0 = Narrowphase.circleHeightfield_Ray.v0,
      v1 = Narrowphase.circleHeightfield_v1;
    // Get the index of the points to test against
    let idxA = Math.floor((circlePos[0] - radius - hfPos[0]) / w),
      idxB = Math.ceil((circlePos[0] + radius - hfPos[0]) / w);
    /*if(idxB < 0 || idxA >= data.length)
        return justTest ? false : 0;*/
    if (idxA < 0) {
      idxA = 0;
    }
    if (idxB >= data.length) {
      idxB = data.length - 1;
    }
    // Get max and min
    let max = data[idxA],
      min = data[idxB];
    for (let i = idxA; i < idxB; i++) {
      if (data[i] < min) {
        min = data[i];
      }
      if (data[i] > max) {
        max = data[i];
      }
    }
    if (circlePos[1] - radius > max) {
      return justTest ? false : 0;
    }
    /*
    if(circlePos[1]+radius < min){
        // Below the minimum point... We can just guess.
        // TODO
    }
    */
    // 1. Check so center of circle is not inside the field. If it is, this wont work...
    // 2. For each edge
    // 2. 1. Get point on circle that is closest to the edge (scale normal with -radius)
    // 2. 2. Check if point is inside.
    let found = false;
    // Check all edges first
    for (let i = idxA; i < idxB; i++) {
      // Get points
      vec2.set(Ray.v0, i * w, data[i]);
      vec2.set(v1, (i + 1) * w, data[i + 1]);
      vec2.add(Ray.v0, Ray.v0, hfPos);
      vec2.add(v1, v1, hfPos);
      // Get normal
      vec2.sub(worldNormal, v1, Ray.v0);
      vec2.rotate(worldNormal, worldNormal, Math.PI / 2);
      vec2.normalize(worldNormal, worldNormal);
      // Get point on circle, closest to the edge
      vec2.scale(candidate, worldNormal, -radius);
      vec2.add(candidate, candidate, circlePos);
      // Distance from Ray.v0 to the candidate point
      vec2.sub(dist, candidate, Ray.v0);
      // Check if it is in the element "stick"
      let d = vec2.dot(dist, worldNormal);
      if (candidate[0] >= Ray.v0[0] && candidate[0] < v1[0] && d <= 0) {
        if (justTest) {
          return true;
        }
        found = true;
        // Store the candidate point, projected to the edge
        vec2.scale(dist, worldNormal, -d);
        vec2.add(minCandidate, candidate, dist);
        vec2.copy(minCandidateNormal, worldNormal);
        let c = this.createContactEquation(hfBody, circleBody, hfShape, circleShape);
        // Normal is out of the heightfield
        vec2.copy(c.normalA, minCandidateNormal);
        // Vector from circle to heightfield
        vec2.scale(c.contactPointB, c.normalA, -radius);
        add(c.contactPointB, c.contactPointB, circlePos);
        sub(c.contactPointB, c.contactPointB, circleBody.position);
        vec2.copy(c.contactPointA, minCandidate);
        vec2.sub(c.contactPointA, c.contactPointA, hfBody.position);
        this.contactEquations.push(c);
        if (this.enableFriction) {
          this.frictionEquations.push(this.createFrictionFromContact(c));
        }
      }
    }
    // Check all vertices
    found = false;
    if (radius > 0) {
      for (let i = idxA; i <= idxB; i++) {
        // Get point
        vec2.set(Ray.v0, i * w, data[i]);
        vec2.add(Ray.v0, Ray.v0, hfPos);
        vec2.sub(dist, circlePos, Ray.v0);
        if (vec2.squaredLength(dist) < Math.pow(radius, 2)) {
          if (justTest) {
            return true;
          }
          found = true;
          let c = this.createContactEquation(hfBody, circleBody, hfShape, circleShape);
          // Construct normal - out of heightfield
          vec2.copy(c.normalA, dist);
          vec2.normalize(c.normalA, c.normalA);
          vec2.scale(c.contactPointB, c.normalA, -radius);
          add(c.contactPointB, c.contactPointB, circlePos);
          sub(c.contactPointB, c.contactPointB, circleBody.position);
          sub(c.contactPointA, Ray.v0, hfPos);
          add(c.contactPointA, c.contactPointA, hfPos);
          sub(c.contactPointA, c.contactPointA, hfBody.position);
          this.contactEquations.push(c);
          if (this.enableFriction) {
            this.frictionEquations.push(this.createFrictionFromContact(c));
          }
        }
      }
    }
    if (found) {
      return 1;
    }
    return 0;
  }
  static convexHeightfield_Ray.v0 = vec2.create();
  static convexHeightfield_v1 = vec2.create();
  static convexHeightfield_tilePos = vec2.create();
  static convexHeightfield_tempConvexShape = new Convex({ vertices: [vec2.create(), vec2.create(), vec2.create(), vec2.create()] });
  /**@param {Body}           bi
  * @param {Circle}         si
  * @param {Array}          xi
  * @param {Body}           bj
  * @param {Heightfield}    sj
  * @param {Array}          xj
  * @param {Number}         aj
  */
  convexHeightfield(convexBody, convexShape, convexPos, convexAngle,
    hfBody, hfShape, hfPos, hfAngle, justTest) {
    let data = hfShape.heights,
      w = hfShape.elementWidth,
      Ray.v0 = Narrowphase.convexHeightfield_Ray.v0,
      v1 = Narrowphase.convexHeightfield_v1,
      tilePos = Narrowphase.convexHeightfield_tilePos,
      tileConvex = Narrowphase.convexHeightfield_tempConvexShape;
    // Get the index of the points to test against
    let idxA = Math.floor((convexBody.aabb.lowerBound[0] - hfPos[0]) / w),
      idxB = Math.ceil((convexBody.aabb.upperBound[0] - hfPos[0]) / w);
    if (idxA < 0) {
      idxA = 0;
    }
    if (idxB >= data.length) {
      idxB = data.length - 1;
    }
    // Get max and min
    let max = data[idxA],
      min = data[idxB];
    for (let i = idxA; i < idxB; i++) {
      if (data[i] < min) {
        min = data[i];
      }
      if (data[i] > max) {
        max = data[i];
      }
    }
    if (convexBody.aabb.lowerBound[1] > max) {
      return justTest ? false : 0;
    }
    let found = false;
    let numContacts = 0;
    // Loop over all edges
    // TODO: If possible, construct a convex from several data points (need o check if the points make a convex shape)
    for (let i = idxA; i < idxB; i++) {
      // Get points
      vec2.set(Ray.v0, i * w, data[i]);
      vec2.set(v1, (i + 1) * w, data[i + 1]);
      vec2.add(Ray.v0, Ray.v0, hfPos);
      vec2.add(v1, v1, hfPos);
      // Construct a convex
      let tileHeight = 100; // todo
      vec2.set(tilePos, (v1[0] + Ray.v0[0]) * 0.5, (v1[1] + Ray.v0[1] - tileHeight) * 0.5);
      vec2.sub(tileConvex.vertices[0], v1, tilePos);
      vec2.sub(tileConvex.vertices[1], Ray.v0, tilePos);
      vec2.copy(tileConvex.vertices[2], tileConvex.vertices[1]);
      vec2.copy(tileConvex.vertices[3], tileConvex.vertices[0]);
      tileConvex.vertices[2][1] -= tileHeight;
      tileConvex.vertices[3][1] -= tileHeight;
      // Do convex collision
      numContacts += this.convexConvex(convexBody, convexShape, convexPos, convexAngle,
        hfBody, tileConvex, tilePos, 0, justTest);
    }
    return numContacts;
  }
}
export class Ray {
  /**A line with a start and end point that is used to intersect shapes. For an example, see {{#crossLink "World/raycast:method"}}World.raycast{{/crossLink}}
  * @param {{from, to, checkCollisionResponse, skipBackfaces, collisionMask, collisionGroup, mode, callback}} options
  */
  constructor(options = {}) {
    this.from = options.from ? vec2.fromValues(options.from[0], options.from[1]) : vec2.create();
    this.to = options.to ? vec2.fromValues(options.to[0], options.to[1]) : vec2.create();
    this.checkCollisionResponse = options.checkCollisionResponse !== undefined ? options.checkCollisionResponse : true;
    this.skipBackfaces = !!options.skipBackfaces;
    this.collisionMask = options.collisionMask !== undefined ? options.collisionMask : -1;
    this.collisionGroup = options.collisionGroup !== undefined ? options.collisionGroup : -1;
    this.mode = options.mode !== undefined ? options.mode : Ray.ANY;
    this.callback = options.callback || function (result) { }
    this.direction = vec2.create();
    this.length = 1;
    this.update();
  }
  static CLOSEST = 1;
  static ANY = 2;
  static ALL = 4;
  /**Should be called if you change the from or to point.
  */
  update() {
    // Update .direction and .length
    let d = this.direction;
    vec2.sub(d, this.to, this.from);
    this.length = vec2.length(d);
    vec2.normalize(d, d);
  }
  /**@param {Array} bodies An array of Body objects.
  */
  intersectBodies(result, bodies) {
    for (let i = 0, l = bodies.length; !result.shouldStop(this) && i < l; i++) {
      let body = bodies[i];
      let aabb = body.getAABB();
      if (aabb.overlapsRay(this) >= 0 || aabb.containsPoint(this.from)) {
        this.intersectBody(result, body);
      }
    }
  }
  static intersectBody_worldPosition = vec2.create();
  /**Shoot a ray at a body, get back information about the hit
  * @param {Body} body
  */
  intersectBody(result, body) {
    let checkCollisionResponse = this.checkCollisionResponse;
    if (checkCollisionResponse && !body.collisionResponse) {
      return;
    }
    let worldPosition = Ray.intersectBody_worldPosition;
    for (let i = 0, N = body.shapes.length; i < N; i++) {
      let shape = body.shapes[i];
      if (checkCollisionResponse && !shape.collisionResponse) {
        continue; // Skip
      }
      if ((this.collisionGroup & shape.collisionMask) === 0 || (shape.collisionGroup & this.collisionMask) === 0) {
        continue;
      }
      // Get world angle and position of the shape
      vec2.rotate(worldPosition, shape.position, body.angle);
      vec2.add(worldPosition, worldPosition, body.position);
      let worldAngle = shape.angle + body.angle;
      this.intersectShape(
        result,
        shape,
        worldAngle,
        worldPosition,
        body
      );
      if (result.shouldStop(this)) {
        break;
      }
    }
  }
  /**
  * @param {Shape} shape
  * @param {number} angle
  * @param {array} position
  * @param {Body} body
  */
  intersectShape(result, shape, angle, position, body) {
    let from = this.from;
    // Checking radius
    let distance = Ray.distanceFromIntersectionSquared(from, this.direction, position);
    if (distance > shape.boundingRadius * shape.boundingRadius) {
      return;
    }
    this._currentBody = body;
    this._currentShape = shape;
    shape.raycast(result, this, position, angle);
    this._currentBody = this._currentShape = null;
  }
  /**Get the AABB of the ray.
  * @param {AABB} aabb
  */
  getAABB(result) {
    let to = this.to;
    let from = this.from;
    vec2.set(
      result.lowerBound,
      Math.min(to[0], from[0]),
      Math.min(to[1], from[1])
    );
    vec2.set(
      result.upperBound,
      Math.max(to[0], from[0]),
      Math.max(to[1], from[1])
    );
  }
  static hitPointWorld = vec2.create();
  /**@param {number} fraction
    * @param {array} normal
    * @param {number} [faceIndex=-1]
    * @returns {Boolean} True if the intersections should continue
  */
  reportIntersection(result, fraction, normal, faceIndex) {
    let from = this.from;
    let to = this.to;
    let shape = this._currentShape;
    let body = this._currentBody;
    // Skip back faces?
    if (this.skipBackfaces && vec2.dot(normal, this.direction) > 0) {
      return;
    }
    switch (this.mode) {
      case Ray.ALL:
        result.set(
          normal,
          shape,
          body,
          fraction,
          faceIndex
        );
        this.callback(result);
        break;
      case Ray.CLOSEST:
        // Store if closer than current closest
        if (fraction < result.fraction || !result.hasHit()) {
          result.set(
            normal,
            shape,
            body,
            fraction,
            faceIndex
          );
        }
        break;
      case Ray.ANY:
        // Report and stop.
        result.set(
          normal,
          shape,
          body,
          fraction,
          faceIndex
        );
        break;
    }
  }
  static v0 = vec2.create();
  static intersect = vec2.create();
  static distanceFromIntersectionSquared(from, direction, position) {
    // Ray.v0 is vector from from to position
    vec2.sub(Ray.v0, position, from);
    let dot = vec2.dot(Ray.v0, direction);
    // intersect = direction * dot + from
    vec2.scale(Ray.intersect, direction, dot);
    vec2.add(Ray.intersect, Ray.intersect, from);
    return vec2.squaredDistance(position, Ray.intersect);
  }
}

export class RaycastResult {
  constructor() {
    this.normal = vec2.create();
    this.shape = null;
    this.body = null;
    this.faceIndex = -1;
    this.fraction = -1;
    this.isStopped = false;
  }
  /**Reset all result data. Must be done before re-using the result object.
  */
  reset() {
    vec2.set(this.normal, 0, 0);
    this.shape = null;
    this.body = null;
    this.faceIndex = -1;
    this.fraction = -1;
    this.isStopped = false;
  }
  /**Get the distance to the hit point.
  * @param {Ray} ray
  */
  getHitDistance(ray) {
    return vec2.distance(ray.from, ray.to) * this.fraction;
  }
  /**Returns true if the ray hit something since the last reset().
  */
  hasHit() {
    return this.fraction !== -1;
  }
  /**Get world hit point.
  * @param {array} out
  * @param {Ray} ray
  */
  getHitPoint(out, ray) {
    vec2.lerp(out, ray.from, ray.to, this.fraction);
  }
  /**Can be called while iterating over hits to stop searching for hit points.
  */
  stop() {
    this.isStopped = true;
  }
  /**@private
  * @param {Ray} ray
  * @returns {Boolean}
  */
  shouldStop(ray) {
    return this.isStopped || (this.fraction !== -1 && ray.mode === Ray.ANY);
  }
  /**@private
  * @param {array} normal
  * @param {Shape} shape
  * @param {Body} body
  * @param {number} fraction
  */
  set(
    normal,
    shape,
    body,
    fraction,
    faceIndex
  ) {
    vec2.copy(this.normal, normal);
    this.shape = shape;
    this.body = body;
    this.fraction = fraction;
    this.faceIndex = faceIndex;
  }
}

export class SAPBroadphase extends Broadphase {
  constructor() {
    super(Broadphase.SAP);
    /**List of bodies currently in the broadphase.*/
    this.axisList = [];
    /**The axis to sort along. 0 means x-axis and 1 y-axis. If your bodies are more spread out over the X axis, set axisIndex to 0, and you will gain some performance*/
    this.axisIndex = 0;
    let that = this;
    this._addBodyHandler = (e) => {
      that.axisList.push(e.body);
    }
    this._removeBodyHandler = (e) => {
      // Remove from list
      let idx = that.axisList.indexOf(e.body);
      if (idx !== -1) {
        that.axisList.splice(idx, 1);
      }
    }
  }
  /**Change the world
  * @param {World} world
  */
  setWorld(world) {
    // Clear the old axis array
    this.axisList.length = 0;
    // Add all bodies from the new world
    Utils.appendArray(this.axisList, world.bodies);
    // Remove old handlers, if any
    world
      .off("addBody", this._addBodyHandler)
      .off("removeBody", this._removeBodyHandler);
    // Add handlers to update the list of bodies.
    world.on("addBody", this._addBodyHandler).on("removeBody", this._removeBodyHandler);
    this.world = world;
  }
  /**Sorts bodies along an axis.
  * @param {Array} a
  * @param {number} axisIndex
  * @returns {Array}
  */
  static sortAxisList(a, axisIndex) {
    axisIndex = axisIndex | 0;
    for (let i = 1, l = a.length; i < l; i++) {
      let v = a[i];
      for (let j = i - 1; j >= 0; j--) {
        if (a[j].aabb.lowerBound[axisIndex] <= v.aabb.lowerBound[axisIndex]) {
          break;
        }
        a[j + 1] = a[j];
      }
      a[j + 1] = v;
    }
    return a;
  }
  sortList() {
    let bodies = this.axisList,
      axisIndex = this.axisIndex;
    // Sort the lists
    SAPBroadphase.SAPBroadphase.sortAxisList(bodies, axisIndex);
  }
  /**Get the colliding pairs
  * @param {World} world
  * @returns {Array}
  */
  getCollisionPairs(world) {
    let bodies = this.axisList,
      result = this.result,
      axisIndex = this.axisIndex;
    result.length = 0;
    // Update all AABBs if needed
    let l = bodies.length;
    while (l--) {
      let b = bodies[l];
      if (b.aabbNeedsUpdate) {
        b.updateAABB();
      }
    }
    // Sort the lists
    this.sortList();
    // Look through the X list
    for (let i = 0, N = bodies.length | 0; i !== N; i++) {
      let bi = bodies[i];
      for (let j = i + 1; j < N; j++) {
        let bj = bodies[j];
        // Bounds overlap?
        let overlaps = (bj.aabb.lowerBound[axisIndex] <= bi.aabb.upperBound[axisIndex]);
        if (!overlaps) {
          break;
        }
        if (Broadphase.canCollide(bi, bj) && this.boundingVolumeCheck(bi, bj)) {
          result.push(bi, bj);
        }
      }
    }
    return result;
  }
  /**Returns all the bodies within an AABB.
  * @param {World} world
  * @param {AABB} aabb
  * @param {array} result An array to store resulting bodies in.
  * @returns {array}
  */
  aabbQuery(world, aabb, result) {
    result = result || [];
    this.sortList();
    let axisIndex = this.axisIndex;
    let axis = 'x';
    if (axisIndex === 1) { axis = 'y'; }
    if (axisIndex === 2) { axis = 'z'; }
    let axisList = this.axisList;
    let lower = aabb.lowerBound[axis];
    let upper = aabb.upperBound[axis];
    for (let i = 0; i < axisList.length; i++) {
      let b = axisList[i];
      if (b.aabbNeedsUpdate) {
        b.updateAABB();
      }
      if (b.aabb.overlaps(aabb)) {
        result.push(b);
      }
    }
    return result;
  }
}
export class Constraint {
  /**Base constraint class.
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {Number} type
  * @param {Object} [options]
  * @param {Object} [options.collideConnected=true]
  */
  constructor(bodyA, bodyB, type, options) {
    /**The type of constraint. May be one of Constraint.DISTANCE, Constraint.GEAR, Constraint.LOCK, Constraint.PRISMATIC or Constraint.REVOLUTE.
    * @property {number} type
    */
    this.type = type;
    options = Utils.defaults(options, {
      collideConnected: true,
      wakeUpBodies: true,
    });
    this.equations = [];
    this.bodyA = bodyA;
    this.bodyB = bodyB;
    this.collideConnected = options.collideConnected;
    // Wake up bodies when connected
    if (options.wakeUpBodies) {
      if (bodyA) {
        bodyA.wakeUp();
      }
      if (bodyB) {
        bodyB.wakeUp();
      }
    }
  }
  /**Updates the internal constraint parameters before solve.
  */
  update() {
    throw new Error("method update() not implmemented in this Constraint subclass!");
  }
  static DISTANCE = 1;
  static GEAR = 2;
  static LOCK = 3;
  static PRISMATIC = 4;
  static REVOLUTE = 5;
  /**Set stiffness for this constraint.
  * @param {Number} stiffness
  */
  setStiffness(stiffness) {
    let eqs = this.equations;
    for (let i = 0; i !== eqs.length; i++) {
      let eq = eqs[i];
      eq.stiffness = stiffness;
      eq.needsUpdate = true;
    }
  }
  /**Set relaxation for this constraint.
  * @param {Number} relaxation
  */
  setRelaxation(relaxation) {
    let eqs = this.equations;
    for (let i = 0; i !== eqs.length; i++) {
      let eq = eqs[i];
      eq.relaxation = relaxation;
      eq.needsUpdate = true;
    }
  }
}

export class DistanceConstraint extends Constraint {
  /**Constraint that tries to keep the distance between two bodies constant
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {{distance, localAnchorA, localAnchorB, maxForce}} options
  * @example
  * // If distance is not given as an option, then the current distance between the bodies is used.
  * // In this example, the bodies will be constrained to have a distance of 2 between their centers.
  * let bodyA = new Body({ mass: 1, position: [-1, 0] });
  * let bodyB = new Body({ mass: 1, position: [1, 0] });
  * let constraint = new DistanceConstraint(bodyA, bodyB);
  * world.addConstraint(constraint);
  * @example
  * // Manually set the distance and anchors
  * let constraint = new DistanceConstraint(bodyA, bodyB, {
  *   distance: 1,          // Distance to keep between the points
  *   localAnchorA: [1, 0], // Point on bodyA
  *   localAnchorB: [-1, 0] // Point on bodyB
  * });
  * world.addConstraint(constraint);
  */
  constructor(bodyA, bodyB, options = Utils.defaults(options, { localAnchorA: [0, 0], localAnchorB: [0, 0] })) {
    super(bodyA, bodyB, Constraint.DISTANCE, options);
    this.localAnchorA = vec2.fromValues(options.localAnchorA[0], options.localAnchorA[1]);
    this.localAnchorB = vec2.fromValues(options.localAnchorB[0], options.localAnchorB[1]);
    let localAnchorA = this.localAnchorA;
    let localAnchorB = this.localAnchorB;
    this.distance = 0;
    if (typeof (options.distance) === 'number') {
      this.distance = options.distance;
    } else {
      // Use the current world distance between the world anchor points.
      let worldAnchorA = vec2.create(),
        worldAnchorB = vec2.create(),
        r = vec2.create();
      // Transform local anchors to world
      vec2.rotate(worldAnchorA, localAnchorA, bodyA.angle);
      vec2.rotate(worldAnchorB, localAnchorB, bodyB.angle);
      vec2.add(r, bodyB.position, worldAnchorB);
      vec2.sub(r, r, worldAnchorA);
      vec2.sub(r, r, bodyA.position);
      this.distance = vec2.length(r);
    }
    let maxForce;
    if (typeof (options.maxForce) === "undefined") {
      maxForce = Number.MAX_VALUE;
    } else {
      maxForce = options.maxForce;
    }
    let normal = new Equation(bodyA, bodyB, -maxForce, maxForce); // Just in the normal direction
    this.equations = [normal];
    this.maxForce = maxForce;
    let r = vec2.create();
    let ri = vec2.create(); // worldAnchorA
    let rj = vec2.create(); // worldAnchorB
    let that = this;
    normal.computeGq = function () {
      let bodyA = this.bodyA,
        bodyB = this.bodyB,
        xi = bodyA.position,
        xj = bodyB.position;
      // Transform local anchors to world
      vec2.rotate(ri, localAnchorA, bodyA.angle);
      vec2.rotate(rj, localAnchorB, bodyB.angle);
      vec2.add(r, xj, rj);
      vec2.sub(r, r, ri);
      vec2.sub(r, r, xi);
      //vec2.sub(r, bodyB.position, bodyA.position);
      return vec2.length(r) - that.distance;
    }
    // Make the contact constraint bilateral
    this.setMaxForce(maxForce);
    this.upperLimitEnabled = false;
    this.upperLimit = 1;
    this.lowerLimitEnabled = false;
    this.lowerLimit = 0;
    this.position = 0;
  }
  static n = vec2.create();
  static ri = vec2.create(); // worldAnchorA
  static rj = vec2.create(); // worldAnchorB
  update() {
    let normal = this.equations[0],
      bodyA = this.bodyA,
      bodyB = this.bodyB,
      distance = this.distance,
      xi = bodyA.position,
      xj = bodyB.position,
      normalEquation = this.equations[0],
      G = normal.G;
    // Transform local anchors to world
    vec2.rotate(ri, this.localAnchorA, bodyA.angle);
    vec2.rotate(rj, this.localAnchorB, bodyB.angle);
    // Get world anchor points and normal
    vec2.add(n, xj, rj);
    vec2.sub(n, n, ri);
    vec2.sub(n, n, xi);
    this.position = vec2.length(n);
    let violating = false;
    if (this.upperLimitEnabled) {
      if (this.position > this.upperLimit) {
        normalEquation.maxForce = 0;
        normalEquation.minForce = -this.maxForce;
        this.distance = this.upperLimit;
        violating = true;
      }
    }
    if (this.lowerLimitEnabled) {
      if (this.position < this.lowerLimit) {
        normalEquation.maxForce = this.maxForce;
        normalEquation.minForce = 0;
        this.distance = this.lowerLimit;
        violating = true;
      }
    }
    if ((this.lowerLimitEnabled || this.upperLimitEnabled) && !violating) {
      // No constraint needed.
      normalEquation.enabled = false;
      return;
    }
    normalEquation.enabled = true;
    vec2.normalize(n, n);
    // Caluclate cross products
    let rixn = vec2.crossLength(ri, n),
      rjxn = vec2.crossLength(rj, n);
    // G = [-n -rixn n rjxn]
    G[0] = -n[0];
    G[1] = -n[1];
    G[2] = -rixn;
    G[3] = n[0];
    G[4] = n[1];
    G[5] = rjxn;
  }
  /**Set the max force to be used
  * @param {Number} maxForce
  */
  setMaxForce(maxForce) {
    let normal = this.equations[0];
    normal.minForce = -maxForce;
    normal.maxForce = maxForce;
  }
  /**Get the max force
  * @returns {Number}
  */
  getMaxForce() {
    let normal = this.equations[0];
    return normal.maxForce;
  }
}

export class GearConstraint extends Constraint {
  /**Constrains the angle of two bodies to each other to be equal. If a gear ratio is not one, the angle of bodyA must be a multiple of the angle of bodyB.
  * @param {Body}            bodyA
  * @param {Body}            bodyB
  * @param {{angle, ratio, maxTorque}} options
  * @example
  *     let constraint = new GearConstraint(bodyA, bodyB);
  *     world.addConstraint(constraint);
  * @example
  *     let constraint = new GearConstraint(bodyA, bodyB, {
  *         ratio: 2,
  *         maxTorque: 1000
  *     });
  *     world.addConstraint(constraint);
  */
  constructor(bodyA, bodyB, options = { ratio: 1, angle: undefined }) {
    super(bodyA, bodyB, Constraint.GEAR, options);
    this.ratio = options.ratio;
    if (options.angle === undefined) {
      options.angle = bodyB.angle - this.ratio * bodyA.angle;
    }
    this.ratio = options.ratio;
    this.equations = [
      new AngleLockEquation(bodyA, bodyB, options),
    ];
    // Set max torque
    if (options.maxTorque !== undefined) {
      this.setMaxTorque(options.maxTorque);
    }
  }
  update() {
    let eq = this.equations[0];
    if (eq.ratio !== this.ratio) {
      eq.setRatio(this.ratio);
    }
    eq.angle = this.angle;
  }
  /**Set the max torque for the constraint.
  * @param {Number} torque
  */
  setMaxTorque(torque) {
    this.equations[0].setMaxTorque(torque);
  }
  /**Get the max torque for the constraint.
  * @returns {Number}
  */
  getMaxTorque(torque) {
    return this.equations[0].maxForce;
  }
}
export class LockConstraint extends Constraint {
  /**Locks the relative position and rotation between two bodies
  * @class LockConstraint
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {{localOffsetB, localAngleB, maxForce}} options
  * @example
  *     // Locks the relative position and rotation between bodyA and bodyB
  *     let constraint = new LockConstraint(bodyA, bodyB);
  *     world.addConstraint(constraint);
  */
  constructor(bodyA, bodyB, options = { localAngleB=0, maxForce: Number.MAX_VALUE }) {
    super(bodyA, bodyB, Constraint.LOCK, options);
    let maxForce = options.maxForce;
    let localAngleB = options.localAngleB;
    let x = new Equation(bodyA, bodyB, -maxForce, maxForce),
      y = new Equation(bodyA, bodyB, -maxForce, maxForce),
      rot = new Equation(bodyA, bodyB, -maxForce, maxForce);
    let l = vec2.create(),
      g = vec2.create(),
      that = this;
    x.computeGq = function () {
      vec2.rotate(l, that.localOffsetB, bodyA.angle);
      vec2.sub(g, bodyB.position, bodyA.position);
      vec2.sub(g, g, l);
      return g[0];
    }
    y.computeGq = function () {
      vec2.rotate(l, that.localOffsetB, bodyA.angle);
      vec2.sub(g, bodyB.position, bodyA.position);
      vec2.sub(g, g, l);
      return g[1];
    }
    let r = vec2.create(),
      t = vec2.create();
    rot.computeGq = function () {
      vec2.rotate(r, that.localOffsetB, bodyB.angle - that.localAngleB);
      vec2.scale(r, r, -1);
      vec2.sub(g, bodyA.position, bodyB.position);
      vec2.add(g, g, r);
      vec2.rotate(t, r, -Math.PI / 2);
      vec2.normalize(t, t);
      return vec2.dot(g, t);
    }
    this.localOffsetB = vec2.create();
    if (options.localOffsetB) {
      vec2.copy(this.localOffsetB, options.localOffsetB);
    } else {
      // Construct from current positions
      vec2.sub(this.localOffsetB, bodyB.position, bodyA.position);
      vec2.rotate(this.localOffsetB, this.localOffsetB, -bodyA.angle);
    }
    this.localAngleB = 0;
    if (typeof (options.localAngleB) === 'number') {
      this.localAngleB = options.localAngleB;
    } else {
      // Construct
      this.localAngleB = bodyB.angle - bodyA.angle;
    }
    this.equations.push(x, y, rot);
    this.setMaxForce(maxForce);
  }
  /**Set the maximum force to be applied.
  * @param {Number} force
  */
  setMaxForce(force) {
    let eqs = this.equations;
    for (let i = 0; i < this.equations.length; i++) {
      eqs[i].maxForce = force;
      eqs[i].minForce = -force;
    }
  }
  /**Get the max force.
  * @returns {Number}
  */
  getMaxForce() {
    return this.equations[0].maxForce;
  }
  static l = vec2.create();
  static r = vec2.create();
  static t = vec2.create();
  static xAxis = vec2.fromValues(1, 0);
  static yAxis = vec2.fromValues(0, 1);
  update() {
    let x = this.equations[0],
      y = this.equations[1],
      rot = this.equations[2],
      bodyA = this.bodyA,
      bodyB = this.bodyB;
    vec2.rotate(l, this.localOffsetB, bodyA.angle);
    vec2.rotate(r, this.localOffsetB, bodyB.angle - this.localAngleB);
    vec2.scale(r, r, -1);
    vec2.rotate(t, r, Math.PI / 2);
    vec2.normalize(t, t);
    x.G[0] = -1;
    x.G[1] = 0;
    x.G[2] = -vec2.crossLength(l, xAxis);
    x.G[3] = 1;
    y.G[0] = 0;
    y.G[1] = -1;
    y.G[2] = -vec2.crossLength(l, yAxis);
    y.G[4] = 1;
    rot.G[0] = -t[0];
    rot.G[1] = -t[1];
    rot.G[3] = t[0];
    rot.G[4] = t[1];
    rot.G[5] = vec2.crossLength(r, t);
  }
}

export class PrismaticConstraint extends Constraint {
  /**Constraint that only allows bodies to move along a line, relative to each other. See <a href="http://www.iforce2d.net/b2dtut/joints-prismatic">this tutorial</a>. Also called "slider constraint".@param {Body}    bodyA
  * @param {Body}    bodyB
  * @param {{maxForce,localAnchorA,localAnchorB,localAxisA,disableRotationalLock,upperLimit,lowerLimit}} options
  * @todo Ability to create using only a point and a worldAxis
  */
  constructor(bodyA, bodyB, options = {}) {
    super(bodyA, bodyB, Constraint.PRISMATIC, options);
    // Get anchors
    let localAnchorA = vec2.fromValues(0, 0),
      localAxisA = vec2.fromValues(1, 0),
      localAnchorB = vec2.fromValues(0, 0);
    if (options.localAnchorA) { vec2.copy(localAnchorA, options.localAnchorA); }
    if (options.localAxisA) { vec2.copy(localAxisA, options.localAxisA); }
    if (options.localAnchorB) { vec2.copy(localAnchorB, options.localAnchorB); }
    this.localAnchorA = localAnchorA;
    this.localAnchorB = localAnchorB;
    this.localAxisA = localAxisA;
    let maxForce = this.maxForce = typeof (options.maxForce) !== "undefined" ? options.maxForce : Number.MAX_VALUE;
    // Translational part
    let trans = new Equation(bodyA, bodyB, -maxForce, maxForce);
    let ri = new vec2.create(),
      rj = new vec2.create(),
      gg = new vec2.create(),
      t = new vec2.create();
    trans.computeGq = function () {
      // g = ( xj + rj - xi - ri ) * t
      return vec2.dot(gg, t);
    }
    trans.updateJacobian = function () {
      let G = this.G,
        xi = bodyA.position,
        xj = bodyB.position;
      vec2.rotate(ri, localAnchorA, bodyA.angle);
      vec2.rotate(rj, localAnchorB, bodyB.angle);
      vec2.add(gg, xj, rj);
      vec2.sub(gg, gg, xi);
      vec2.sub(gg, gg, ri);
      vec2.rotate(t, localAxisA, bodyA.angle + Math.PI / 2);
      G[0] = -t[0];
      G[1] = -t[1];
      G[2] = -vec2.crossLength(ri, t) + vec2.crossLength(t, gg);
      G[3] = t[0];
      G[4] = t[1];
      G[5] = vec2.crossLength(rj, t);
    }
    this.equations.push(trans);
    // Rotational part
    if (!options.disableRotationalLock) {
      let rot = new RotationalLockEquation(bodyA, bodyB, -maxForce, maxForce);
      this.equations.push(rot);
    }
    this.position = 0;
    // Is this one used at all?
    this.velocity = 0;
    this.lowerLimitEnabled = typeof (options.lowerLimit) !== "undefined" ? true : false;
    this.upperLimitEnabled = typeof (options.upperLimit) !== "undefined" ? true : false;
    this.lowerLimit = typeof (options.lowerLimit) !== "undefined" ? options.lowerLimit : 0;
    this.upperLimit = typeof (options.upperLimit) !== "undefined" ? options.upperLimit : 1;
    // Equations used for limits
    this.upperLimitEquation = new ContactEquation(bodyA, bodyB);
    this.lowerLimitEquation = new ContactEquation(bodyA, bodyB);
    // Set max/min forces
    this.upperLimitEquation.minForce = this.lowerLimitEquation.minForce = 0;
    this.upperLimitEquation.maxForce = this.lowerLimitEquation.maxForce = maxForce;
    this.motorEquation = new Equation(bodyA, bodyB);
    this.motorEnabled = false;
    this.motorSpeed = 0;
    let that = this;
    let motorEquation = this.motorEquation;
    let old = motorEquation.computeGW;
    motorEquation.computeGq = function () { return 0; }
    motorEquation.computeGW = function () {
      let G = this.G,
        bi = this.bodyA,
        bj = this.bodyB,
        vi = bi.velocity,
        vj = bj.velocity,
        wi = bi.angularVelocity,
        wj = bj.angularVelocity;
      return this.gmult(G, vi, wi, vj, wj) + that.motorSpeed;
    }
  }

  static worldAxisA = vec2.create();
  static worldAnchorA = vec2.create();
  static worldAnchorB = vec2.create();
  static orientedAnchorA = vec2.create();
  static orientedAnchorB = vec2.create();
  static tmp = vec2.create();
  /**Update the constraint equations. Should be done if any of the bodies changed position, before solving.
  */
  update() {
    let eqs = this.equations,
      trans = eqs[0],
      upperLimit = this.upperLimit,
      lowerLimit = this.lowerLimit,
      upperLimitEquation = this.upperLimitEquation,
      lowerLimitEquation = this.lowerLimitEquation,
      bodyA = this.bodyA,
      bodyB = this.bodyB,
      localAxisA = this.localAxisA,
      localAnchorA = this.localAnchorA,
      localAnchorB = this.localAnchorB;
    trans.updateJacobian();
    // Transform local things to world
    vec2.rotate(worldAxisA, localAxisA, bodyA.angle);
    vec2.rotate(orientedAnchorA, localAnchorA, bodyA.angle);
    vec2.add(worldAnchorA, orientedAnchorA, bodyA.position);
    vec2.rotate(orientedAnchorB, localAnchorB, bodyB.angle);
    vec2.add(worldAnchorB, orientedAnchorB, bodyB.position);
    let relPosition = this.position = vec2.dot(worldAnchorB, worldAxisA) - vec2.dot(worldAnchorA, worldAxisA);
    // Motor
    if (this.motorEnabled) {
      // G = [ a     a x ri   -a   -a x rj ]
      let G = this.motorEquation.G;
      G[0] = worldAxisA[0];
      G[1] = worldAxisA[1];
      G[2] = vec2.crossLength(worldAxisA, orientedAnchorB);
      G[3] = -worldAxisA[0];
      G[4] = -worldAxisA[1];
      G[5] = -vec2.crossLength(worldAxisA, orientedAnchorA);
    }
    /*
        Limits strategy:
        Add contact equation, with normal along the constraint axis.
        min/maxForce is set so the constraint is repulsive in the correct direction.
        Some offset is added to either equation.contactPointA or .contactPointB to get the correct upper/lower limit.
                 ^
                 |
      upperLimit x
                 |    ------
         anchorB x<---|  B |
                 |    |    |
        ------   |    ------
        |    |   |
        |  A |-->x anchorA
        ------   |
                 x lowerLimit
                 |
                axis
     */

    if (this.upperLimitEnabled && relPosition > upperLimit) {
      // Update contact constraint normal, etc
      vec2.scale(upperLimitEquation.normalA, worldAxisA, -1);
      vec2.sub(upperLimitEquation.contactPointA, worldAnchorA, bodyA.position);
      vec2.sub(upperLimitEquation.contactPointB, worldAnchorB, bodyB.position);
      vec2.scale(tmp, worldAxisA, upperLimit);
      vec2.add(upperLimitEquation.contactPointA, upperLimitEquation.contactPointA, tmp);
      if (eqs.indexOf(upperLimitEquation) === -1) {
        eqs.push(upperLimitEquation);
      }
    } else {
      let idx = eqs.indexOf(upperLimitEquation);
      if (idx !== -1) {
        eqs.splice(idx, 1);
      }
    }
    if (this.lowerLimitEnabled && relPosition < lowerLimit) {
      // Update contact constraint normal, etc
      vec2.scale(lowerLimitEquation.normalA, worldAxisA, 1);
      vec2.sub(lowerLimitEquation.contactPointA, worldAnchorA, bodyA.position);
      vec2.sub(lowerLimitEquation.contactPointB, worldAnchorB, bodyB.position);
      vec2.scale(tmp, worldAxisA, lowerLimit);
      vec2.sub(lowerLimitEquation.contactPointB, lowerLimitEquation.contactPointB, tmp);
      if (eqs.indexOf(lowerLimitEquation) === -1) {
        eqs.push(lowerLimitEquation);
      }
    } else {
      let idx = eqs.indexOf(lowerLimitEquation);
      if (idx !== -1) {
        eqs.splice(idx, 1);
      }
    }
  }
  /**Enable the motor
  */
  enableMotor() {
    if (this.motorEnabled) {
      return;
    }
    this.equations.push(this.motorEquation);
    this.motorEnabled = true;
  }
  /**Disable the rotational motor
  */
  disableMotor() {
    if (!this.motorEnabled) {
      return;
    }
    let i = this.equations.indexOf(this.motorEquation);
    this.equations.splice(i, 1);
    this.motorEnabled = false;
  }
  /**Set the constraint limits.
  * @param {number} lower Lower limit.
  * @param {number} upper Upper limit.
  */
  setLimits(lower, upper) {
    if (typeof (lower) === 'number') {
      this.lowerLimit = lower;
      this.lowerLimitEnabled = true;
    } else {
      this.lowerLimit = lower;
      this.lowerLimitEnabled = false;
    }
    if (typeof (upper) === 'number') {
      this.upperLimit = upper;
      this.upperLimitEnabled = true;
    } else {
      this.upperLimit = upper;
      this.upperLimitEnabled = false;
    }
  }
}

class RevoluteConstraint extends Constraint {
  let worldPivotA = vec2.create(),
  worldPivotB = vec2.create(),
  xAxis = vec2.fromValues(1, 0),
  yAxis = vec2.fromValues(0, 1),
  g = vec2.create();
/**Connects two bodies at given offset points, letting them rotate relative to each other around this point.
* @class RevoluteConstraint
* @constructor
* @author schteppe
* @param {Body}    bodyA
* @param {Body}    bodyB
* @param {Object}  [options]
* @param {Array}   [options.worldPivot] A pivot point given in world coordinates. If specified, localPivotA and localPivotB are automatically computed from this value.
* @param {Array}   [options.localPivotA] The point relative to the center of mass of bodyA which bodyA is constrained to.
* @param {Array}   [options.localPivotB] See localPivotA.
* @param {Number}  [options.maxForce] The maximum force that should be applied to constrain the bodies.
* @extends Constraint
*
* @example
*     // This will create a revolute constraint between two bodies with pivot point in between them.
*     let bodyA = new Body({ mass: 1, position: [-1, 0] });
*     let bodyB = new Body({ mass: 1, position: [1, 0] });
*     let constraint = new RevoluteConstraint(bodyA, bodyB, {
*         worldPivot: [0, 0]
*     });
*     world.addConstraint(constraint);
*
*     // Using body-local pivot points, the constraint could have been constructed like this:
*     let constraint = new RevoluteConstraint(bodyA, bodyB, {
*         localPivotA: [1, 0],
*         localPivotB: [-1, 0]
*     });
*/
function RevoluteConstraint(bodyA, bodyB, options) {
  options = options || {}
  Constraint.call(this, bodyA, bodyB, Constraint.REVOLUTE, options);
  let maxForce = this.maxForce = typeof (options.maxForce) !== "undefined" ? options.maxForce : Number.MAX_VALUE;
  /**@property {Array} pivotA
  */
  this.pivotA = vec2.create();
  /**@property {Array} pivotB
  */
  this.pivotB = vec2.create();
  if (options.worldPivot) {
    // Compute pivotA and pivotB
    vec2.sub(this.pivotA, options.worldPivot, bodyA.position);
    vec2.sub(this.pivotB, options.worldPivot, bodyB.position);
    // Rotate to local coordinate system
    vec2.rotate(this.pivotA, this.pivotA, -bodyA.angle);
    vec2.rotate(this.pivotB, this.pivotB, -bodyB.angle);
  } else {
    // Get pivotA and pivotB
    vec2.copy(this.pivotA, options.localPivotA);
    vec2.copy(this.pivotB, options.localPivotB);
  }
  // Equations to be fed to the solver
  let eqs = this.equations = [
    new Equation(bodyA, bodyB, -maxForce, maxForce),
    new Equation(bodyA, bodyB, -maxForce, maxForce),
  ];
  let x = eqs[0];
  let y = eqs[1];
  let that = this;
  x.computeGq() {
    vec2.rotate(worldPivotA, that.pivotA, bodyA.angle);
    vec2.rotate(worldPivotB, that.pivotB, bodyB.angle);
    vec2.add(g, bodyB.position, worldPivotB);
    vec2.sub(g, g, bodyA.position);
    vec2.sub(g, g, worldPivotA);
    return vec2.dot(g, xAxis);
  }
  y.computeGq() {
    vec2.rotate(worldPivotA, that.pivotA, bodyA.angle);
    vec2.rotate(worldPivotB, that.pivotB, bodyB.angle);
    vec2.add(g, bodyB.position, worldPivotB);
    vec2.sub(g, g, bodyA.position);
    vec2.sub(g, g, worldPivotA);
    return vec2.dot(g, yAxis);
  }
  y.minForce = x.minForce = -maxForce;
  y.maxForce = x.maxForce = maxForce;
  this.motorEquation = new RotationalVelocityEquation(bodyA, bodyB);
  /**Indicates whether the motor is enabled. Use .enableMotor() to enable the constraint motor.
  * @property {Boolean} motorEnabled
  * @readOnly
  */
  this.motorEnabled = false;
  /**The constraint position.
  * @property angle
  * @type {Number}
  * @readOnly
  */
  this.angle = 0;
  /**Set to true to enable lower limit
  * @property lowerLimitEnabled
  * @type {Boolean}
  */
  this.lowerLimitEnabled = false;
  /**Set to true to enable upper limit
  * @property upperLimitEnabled
  * @type {Boolean}
  */
  this.upperLimitEnabled = false;
  /**The lower limit on the constraint angle.
  * @property lowerLimit
  * @type {Boolean}
  */
  this.lowerLimit = 0;
  /**The upper limit on the constraint angle.
  * @property upperLimit
  * @type {Boolean}
  */
  this.upperLimit = 0;
  this.upperLimitEquation = new RotationalLockEquation(bodyA, bodyB);
  this.lowerLimitEquation = new RotationalLockEquation(bodyA, bodyB);
  this.upperLimitEquation.minForce = 0;
  this.lowerLimitEquation.maxForce = 0;
}
RevoluteConstraint.prototype = new Constraint();
Revoluteconstructor = RevoluteConstraint;
/**Set the constraint angle limits.
* @param {number} lower Lower angle limit.
* @param {number} upper Upper angle limit.
*/
RevolutesetLimits(lower, upper) {
  if (typeof (lower) === 'number') {
    this.lowerLimit = lower;
    this.lowerLimitEnabled = true;
  } else {
    this.lowerLimit = lower;
    this.lowerLimitEnabled = false;
  }
  if (typeof (upper) === 'number') {
    this.upperLimit = upper;
    this.upperLimitEnabled = true;
  } else {
    this.upperLimit = upper;
    this.upperLimitEnabled = false;
  }
}
Revoluteupdate() {
  let bodyA = this.bodyA,
    bodyB = this.bodyB,
    pivotA = this.pivotA,
    pivotB = this.pivotB,
    eqs = this.equations,
    normal = eqs[0],
    tangent = eqs[1],
    x = eqs[0],
    y = eqs[1],
    upperLimit = this.upperLimit,
    lowerLimit = this.lowerLimit,
    upperLimitEquation = this.upperLimitEquation,
    lowerLimitEquation = this.lowerLimitEquation;
  let relAngle = this.angle = bodyB.angle - bodyA.angle;
  if (this.upperLimitEnabled && relAngle > upperLimit) {
    upperLimitEquation.angle = upperLimit;
    if (eqs.indexOf(upperLimitEquation) === -1) {
      eqs.push(upperLimitEquation);
    }
  } else {
    let idx = eqs.indexOf(upperLimitEquation);
    if (idx !== -1) {
      eqs.splice(idx, 1);
    }
  }
  if (this.lowerLimitEnabled && relAngle < lowerLimit) {
    lowerLimitEquation.angle = lowerLimit;
    if (eqs.indexOf(lowerLimitEquation) === -1) {
      eqs.push(lowerLimitEquation);
    }
  } else {
    let idx = eqs.indexOf(lowerLimitEquation);
    if (idx !== -1) {
      eqs.splice(idx, 1);
    }
  }
  /*
  The constraint violation is
      g = xj + rj - xi - ri
  ...where xi and xj are the body positions and ri and rj world-oriented offset vectors. Differentiate:
      gdot = vj + wj x rj - vi - wi x ri
  We split this into x and y directions. (let x and y be unit vectors along the respective axes)
      gdot * x = ( vj + wj x rj - vi - wi x ri ) * x
               = ( vj*x + (wj x rj)*x -vi*x -(wi x ri)*x
               = ( vj*x + (rj x x)*wj -vi*x -(ri x x)*wi
               = [ -x   -(ri x x)   x   (rj x x)] * [vi wi vj wj]
               = G*W
  ...and similar for y. We have then identified the jacobian entries for x and y directions:
      Gx = [ x   (rj x x)   -x   -(ri x x)]
      Gy = [ y   (rj x y)   -y   -(ri x y)]
   */
  vec2.rotate(worldPivotA, pivotA, bodyA.angle);
  vec2.rotate(worldPivotB, pivotB, bodyB.angle);
  // todo: these are a bit sparse. We could save some computations on making custom eq.computeGW functions, etc
  x.G[0] = -1;
  x.G[1] = 0;
  x.G[2] = -vec2.crossLength(worldPivotA, xAxis);
  x.G[3] = 1;
  x.G[4] = 0;
  x.G[5] = vec2.crossLength(worldPivotB, xAxis);
  y.G[0] = 0;
  y.G[1] = -1;
  y.G[2] = -vec2.crossLength(worldPivotA, yAxis);
  y.G[3] = 0;
  y.G[4] = 1;
  y.G[5] = vec2.crossLength(worldPivotB, yAxis);
}
/**Enable the rotational motor
*/
RevoluteenableMotor() {
  if (this.motorEnabled) {
    return;
  }
  this.equations.push(this.motorEquation);
  this.motorEnabled = true;
}
/**Disable the rotational motor
*/
RevolutedisableMotor() {
  if (!this.motorEnabled) {
    return;
  }
  let i = this.equations.indexOf(this.motorEquation);
  this.equations.splice(i, 1);
  this.motorEnabled = false;
}
/**Check if the motor is enabled.
* @deprecated use property motorEnabled instead.
* @returns {Boolean}
*/
RevolutemotorIsEnabled() {
  return !!this.motorEnabled;
}
/**Set the speed of the rotational constraint motor
* @param {Number} speed
*/
RevolutesetMotorSpeed(speed) {
  if (!this.motorEnabled) {
    return;
  }
  let i = this.equations.indexOf(this.motorEquation);
  this.equations[i].relativeVelocity = speed;
}
/**Get the speed of the rotational constraint motor
* @returns {Number} The current speed, or false if the motor is not enabled.
*/
RevolutegetMotorSpeed() {
  if (!this.motorEnabled) {
    return false;
  }
  return this.motorEquation.relativeVelocity;
}
}, { "../equations/Equation": 22, "../equations/RotationalLockEquation": 24, "../equations/RotationalVelocityEquation": 25, "../math/vec2": 30, "./Constraint": 14 }], 20: [function (_dereq_, module, exports) {
  let Equation = _dereq_("./Equation"),
    vec2 = _dereq_('../math/vec2');
  module.exports = AngleLockEquation;
  /**Locks the relative angle between two bodies. The constraint tries to keep the dot product between two vectors, local in each body, to zero. The local angle in body i is a parameter.
  *
  * @class AngleLockEquation
  * @constructor
  * @extends Equation
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {Object} [options]
  * @param {Number} [options.angle] Angle to add to the local vector in body A.
  * @param {Number} [options.ratio] Gear ratio
  */
  function AngleLockEquation(bodyA, bodyB, options) {
    options = options || {}
    Equation.call(this, bodyA, bodyB, -Number.MAX_VALUE, Number.MAX_VALUE);
    this.angle = options.angle || 0;
    /**The gear ratio.
    * @property {Number} ratio
    * @private
    * @see setRatio
    */
    this.ratio = typeof (options.ratio) === "number" ? options.ratio : 1;
    this.setRatio(this.ratio);
  }
  AngleLockEquation.prototype = new Equation();
  AngleLockEquation.prototype.constructor = AngleLockEquation;
  AngleLockEquation.prototype.computeGq() {
    return this.ratio * this.bodyA.angle - this.bodyB.angle + this.angle;
  }
  /**Set the gear ratio for this equation
  * @param {Number} ratio
  */
  AngleLockEquation.prototype.setRatio(ratio) {
    let G = this.G;
    G[2] = ratio;
    G[5] = -1;
    this.ratio = ratio;
  }
  /**Set the max force for the equation.
  * @param {Number} torque
  */
  AngleLockEquation.prototype.setMaxTorque(torque) {
    this.maxForce = torque;
    this.minForce = -torque;
  }
}, { "../math/vec2": 30, "./Equation": 22 }], 21: [function (_dereq_, module, exports) {
  let Equation = _dereq_("./Equation"),
    vec2 = _dereq_('../math/vec2');
  module.exports = ContactEquation;
  /**Non-penetration constraint equation. Tries to make the contactPointA and contactPointB vectors coincide, while keeping the applied force repulsive.
  *
  * @class ContactEquation
  * @constructor
  * @extends Equation
  * @param {Body} bodyA
  * @param {Body} bodyB
  */
  function ContactEquation(bodyA, bodyB) {
    Equation.call(this, bodyA, bodyB, 0, Number.MAX_VALUE);
    /**Vector from body i center of mass to the contact point.
    * @property contactPointA
    * @type {Array}
    */
    this.contactPointA = vec2.create();
    this.penetrationVec = vec2.create();
    /**World-oriented vector from body A center of mass to the contact point.
    * @property contactPointB
    * @type {Array}
    */
    this.contactPointB = vec2.create();
    /**The normal vector, pointing out of body i
    * @property normalA
    * @type {Array}
    */
    this.normalA = vec2.create();
    /**The restitution to use (0=no bounciness, 1=max bounciness).
    * @property restitution
    * @type {Number}
    */
    this.restitution = 0;
    /**This property is set to true if this is the first impact between the bodies (not persistant contact).
    * @property firstImpact
    * @type {Boolean}
    * @readOnly
    */
    this.firstImpact = false;
    /**The shape in body i that triggered this contact.
    * @property shapeA
    * @type {Shape}
    */
    this.shapeA = null;
    /**The shape in body j that triggered this contact.
    * @property shapeB
    * @type {Shape}
    */
    this.shapeB = null;
  }
  ContactEquation.prototype = new Equation();
  ContactEquation.prototype.constructor = ContactEquation;
  ContactEquation.prototype.computeB(a, b, h) {
    let bi = this.bodyA,
      bj = this.bodyB,
      ri = this.contactPointA,
      rj = this.contactPointB,
      xi = bi.position,
      xj = bj.position;
    let penetrationVec = this.penetrationVec,
      n = this.normalA,
      G = this.G;
    // Caluclate cross products
    let rixn = vec2.crossLength(ri, n),
      rjxn = vec2.crossLength(rj, n);
    // G = [-n -rixn n rjxn]
    G[0] = -n[0];
    G[1] = -n[1];
    G[2] = -rixn;
    G[3] = n[0];
    G[4] = n[1];
    G[5] = rjxn;
    // Calculate q = xj+rj -(xi+ri) i.e. the penetration vector
    vec2.add(penetrationVec, xj, rj);
    vec2.sub(penetrationVec, penetrationVec, xi);
    vec2.sub(penetrationVec, penetrationVec, ri);
    // Compute iteration
    let GW, Gq;
    if (this.firstImpact && this.restitution !== 0) {
      Gq = 0;
      GW = (1 / b) * (1 + this.restitution) * this.computeGW();
    } else {
      Gq = vec2.dot(n, penetrationVec) + this.offset;
      GW = this.computeGW();
    }
    let GiMf = this.computeGiMf();
    let B = - Gq * a - GW * b - h * GiMf;
    return B;
  }
  let vi = vec2.create();
  let vj = vec2.create();
  let relVel = vec2.create();
  /**Get the relative velocity along the normal vector.
  * @returns {number}
  */
  ContactEquation.prototype.getVelocityAlongNormal() {
    this.bodyA.getVelocityAtPoint(vi, this.contactPointA);
    this.bodyB.getVelocityAtPoint(vj, this.contactPointB);
    vec2.subtract(relVel, vi, vj);
    return vec2.dot(this.normalA, relVel);
  }
}, { "../math/vec2": 30, "./Equation": 22 }], 22: [function (_dereq_, module, exports) {
  module.exports = Equation;
  let vec2 = _dereq_('../math/vec2'),
    Utils = _dereq_('../utils/Utils'),
    Body = _dereq_('../objects/Body');
  /**Base class for constraint equations.
  * @class Equation
  * @constructor
  * @param {Body} bodyA First body participating in the equation
  * @param {Body} bodyB Second body participating in the equation
  * @param {number} minForce Minimum force to apply. Default: -Number.MAX_VALUE
  * @param {number} maxForce Maximum force to apply. Default: Number.MAX_VALUE
  */
  function Equation(bodyA, bodyB, minForce, maxForce) {
    /**Minimum force to apply when solving.
    * @property minForce
    * @type {Number}
    */
    this.minForce = typeof (minForce) === "undefined" ? -Number.MAX_VALUE : minForce;
    /**Max force to apply when solving.
    * @property maxForce
    * @type {Number}
    */
    this.maxForce = typeof (maxForce) === "undefined" ? Number.MAX_VALUE : maxForce;
    /**First body participating in the constraint
    * @property bodyA
    * @type {Body}
    */
    this.bodyA = bodyA;
    /**Second body participating in the constraint
    * @property bodyB
    * @type {Body}
    */
    this.bodyB = bodyB;
    /**The stiffness of this equation. Typically chosen to a large number (~1e7), but can be chosen somewhat freely to get a stable simulation.
    * @property stiffness
    * @type {Number}
    */
    this.stiffness = Equation.DEFAULT_STIFFNESS;
    /**The number of time steps needed to stabilize the constraint equation. Typically between 3 and 5 time steps.
    * @property relaxation
    * @type {Number}
    */
    this.relaxation = Equation.DEFAULT_RELAXATION;
    /**The Jacobian entry of this equation. 6 numbers, 3 per body (x,y,angle).
    * @property G
    * @type {Array}
    */
    this.G = new Utils.ARRAY_TYPE(6);
    for (let i = 0; i < 6; i++) {
      this.G[i] = 0;
    }
    this.offset = 0;
    this.a = 0;
    this.b = 0;
    this.epsilon = 0;
    this.timeStep = 1 / 60;
    /**Indicates if stiffness or relaxation was changed.
    * @property {Boolean} needsUpdate
    */
    this.needsUpdate = true;
    /**The resulting constraint multiplier from the last solve. This is mostly equivalent to the force produced by the constraint.
    * @property multiplier
    * @type {Number}
    */
    this.multiplier = 0;
    /**Relative velocity.
    * @property {Number} relativeVelocity
    */
    this.relativeVelocity = 0;
    /**Whether this equation is enabled or not. If true, it will be added to the solver.
    * @property {Boolean} enabled
    */
    this.enabled = true;
  }
  Equation.prototype.constructor = Equation;
  /**The default stiffness when creating a new Equation.
   * @property {Number} DEFAULT_STIFFNESS
  * @default 1e6
  */
  Equation.DEFAULT_STIFFNESS = 1e6;
  /**The default relaxation when creating a new Equation.
   * @property {Number} DEFAULT_RELAXATION
  * @default 4
  */
  Equation.DEFAULT_RELAXATION = 4;
  /**Compute SPOOK parameters .a, .b and .epsilon according to the current parameters. See equations 9, 10 and 11 in the <a href="http://www8.cs.umu.se/kurser/5DRay.v058/VT09/lectures/spooknotes.pdf">SPOOK notes</a>.
  */
  Equation.prototype.update() {
    let k = this.stiffness,
      d = this.relaxation,
      h = this.timeStep;
    this.a = 4.0 / (h * (1 + 4 * d));
    this.b = (4.0 * d) / (1 + 4 * d);
    this.epsilon = 4.0 / (h * h * k * (1 + 4 * d));
    this.needsUpdate = false;
  }
  /**Multiply a jacobian entry with corresponding positions or velocities
  * @returns {Number}
  */
  Equation.prototype.gmult(G, vi, wi, vj, wj) {
    return G[0] * vi[0] +
      G[1] * vi[1] +
      G[2] * wi +
      G[3] * vj[0] +
      G[4] * vj[1] +
      G[5] * wj;
  }
  /**Computes the RHS of the SPOOK equation
  * @returns {Number}
  */
  Equation.prototype.computeB(a, b, h) {
    let GW = this.computeGW();
    let Gq = this.computeGq();
    let GiMf = this.computeGiMf();
    return - Gq * a - GW * b - GiMf * h;
  }
  /**Computes G\*q, where q are the generalized body coordinates
  * @returns {Number}
  */
  let qi = vec2.create(),
    qj = vec2.create();
  Equation.prototype.computeGq() {
    let G = this.G,
      bi = this.bodyA,
      bj = this.bodyB,
      xi = bi.position,
      xj = bj.position,
      ai = bi.angle,
      aj = bj.angle;
    return this.gmult(G, qi, ai, qj, aj) + this.offset;
  }
  /**Computes G\*W, where W are the body velocities
  * @returns {Number}
  */
  Equation.prototype.computeGW() {
    let G = this.G,
      bi = this.bodyA,
      bj = this.bodyB,
      vi = bi.velocity,
      vj = bj.velocity,
      wi = bi.angularVelocity,
      wj = bj.angularVelocity;
    return this.gmult(G, vi, wi, vj, wj) + this.relativeVelocity;
  }
  /**Computes G\*Wlambda, where W are the body velocities
  * @returns {Number}
  */
  Equation.prototype.computeGWlambda() {
    let G = this.G,
      bi = this.bodyA,
      bj = this.bodyB,
      vi = bi.vlambda,
      vj = bj.vlambda,
      wi = bi.wlambda,
      wj = bj.wlambda;
    return this.gmult(G, vi, wi, vj, wj);
  }
  /**Computes G\*inv(M)\*f, where M is the mass matrix with diagonal blocks for each body, and f are the forces on the bodies.
  * @returns {Number}
  */
  let iMfi = vec2.create(),
    iMfj = vec2.create();
  Equation.prototype.computeGiMf() {
    let bi = this.bodyA,
      bj = this.bodyB,
      fi = bi.force,
      ti = bi.angularForce,
      fj = bj.force,
      tj = bj.angularForce,
      invMassi = bi.invMassSolve,
      invMassj = bj.invMassSolve,
      invIi = bi.invInertiaSolve,
      invIj = bj.invInertiaSolve,
      G = this.G;
    vec2.scale(iMfi, fi, invMassi);
    vec2.multiply(iMfi, bi.massMultiplier, iMfi);
    vec2.scale(iMfj, fj, invMassj);
    vec2.multiply(iMfj, bj.massMultiplier, iMfj);
    return this.gmult(G, iMfi, ti * invIi, iMfj, tj * invIj);
  }
  /**Computes G\*inv(M)\*G'
  * @returns {Number}
  */
  Equation.prototype.computeGiMGt() {
    let bi = this.bodyA,
      bj = this.bodyB,
      invMassi = bi.invMassSolve,
      invMassj = bj.invMassSolve,
      invIi = bi.invInertiaSolve,
      invIj = bj.invInertiaSolve,
      G = this.G;
    return G[0] * G[0] * invMassi * bi.massMultiplier[0] +
      G[1] * G[1] * invMassi * bi.massMultiplier[1] +
      G[2] * G[2] * invIi +
      G[3] * G[3] * invMassj * bj.massMultiplier[0] +
      G[4] * G[4] * invMassj * bj.massMultiplier[1] +
      G[5] * G[5] * invIj;
  }
  let addToWlambda_temp = vec2.create(),
    addToWlambda_Gi = vec2.create(),
    addToWlambda_Gj = vec2.create(),
    addToWlambda_ri = vec2.create(),
    addToWlambda_rj = vec2.create(),
    addToWlambda_Mdiag = vec2.create();
  /**Add constraint velocity to the bodies.
  * @param {Number} deltalambda
  */
  Equation.prototype.addToWlambda(deltalambda) {
    let bi = this.bodyA,
      bj = this.bodyB,
      temp = addToWlambda_temp,
      Gi = addToWlambda_Gi,
      Gj = addToWlambda_Gj,
      ri = addToWlambda_ri,
      rj = addToWlambda_rj,
      invMassi = bi.invMassSolve,
      invMassj = bj.invMassSolve,
      invIi = bi.invInertiaSolve,
      invIj = bj.invInertiaSolve,
      Mdiag = addToWlambda_Mdiag,
      G = this.G;
    Gi[0] = G[0];
    Gi[1] = G[1];
    Gj[0] = G[3];
    Gj[1] = G[4];
    // Add to linear velocity
    // v_lambda += inv(M) * delta_lamba * G
    vec2.scale(temp, Gi, invMassi * deltalambda);
    vec2.multiply(temp, temp, bi.massMultiplier);
    vec2.add(bi.vlambda, bi.vlambda, temp);
    // This impulse is in the offset frame
    // Also add contribution to angular
    //bi.wlambda -= vec2.crossLength(temp,ri);
    bi.wlambda += invIi * G[2] * deltalambda;

    vec2.scale(temp, Gj, invMassj * deltalambda);
    vec2.multiply(temp, temp, bj.massMultiplier);
    vec2.add(bj.vlambda, bj.vlambda, temp);
    //bj.wlambda -= vec2.crossLength(temp,rj);
    bj.wlambda += invIj * G[5] * deltalambda;
  }
  /**Compute the denominator part of the SPOOK equation: C = G\*inv(M)\*G' + eps
  * @param {Number} eps
  * @returns {Number}
  */
  Equation.prototype.computeInvC(eps) {
    return 1.0 / (this.computeGiMGt() + eps);
  }
}, { "../math/vec2": 30, "../objects/Body": 31, "../utils/Utils": 57 }], 23: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2')
    , Equation = _dereq_('./Equation')
    , Utils = _dereq_('../utils/Utils');
  module.exports = FrictionEquation;
  /**Constrains the slipping in a contact along a tangent
  *
  * @class FrictionEquation
  * @constructor
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {Number} slipForce
  * @extends Equation
  */
  function FrictionEquation(bodyA, bodyB, slipForce) {
    Equation.call(this, bodyA, bodyB, -slipForce, slipForce);
    /**Relative vector from center of body A to the contact point, world oriented.
    * @property contactPointA
    * @type {Array}
    */
    this.contactPointA = vec2.create();
    /**Relative vector from center of body B to the contact point, world oriented.
    * @property contactPointB
    * @type {Array}
    */
    this.contactPointB = vec2.create();
    /**Tangent vector that the friction force will act along. World oriented.
    * @property t
    * @type {Array}
    */
    this.t = vec2.create();
    /**ContactEquations connected to this friction equation. The contact equations can be used to rescale the max force for the friction. If more than one contact equation is given, then the max force can be set to the average.
    * @property contactEquations
    * @type {ContactEquation}
    */
    this.contactEquations = [];
    /**The shape in body i that triggered this friction.
    * @property shapeA
    * @type {Shape}
    * @todo Needed? The shape can be looked up via contactEquation.shapeA...
    */
    this.shapeA = null;
    /**The shape in body j that triggered this friction.
    * @property shapeB
    * @type {Shape}
    * @todo Needed? The shape can be looked up via contactEquation.shapeB...
    */
    this.shapeB = null;
    /**The friction coefficient to use.
    * @property frictionCoefficient
    * @type {Number}
    */
    this.frictionCoefficient = 0.3;
  }
  FrictionEquation.prototype = new Equation();
  FrictionEquation.prototype.constructor = FrictionEquation;
  /**Set the slipping condition for the constraint. The friction force cannot be
  * larger than this value.
  * @param {Number} slipForce
  */
  FrictionEquation.prototype.setSlipForce(slipForce) {
    this.maxForce = slipForce;
    this.minForce = -slipForce;
  }
  /**Get the max force for the constraint.
  * @returns {Number}
  */
  FrictionEquation.prototype.getSlipForce() {
    return this.maxForce;
  }
  FrictionEquation.prototype.computeB(a, b, h) {
    let bi = this.bodyA,
      bj = this.bodyB,
      ri = this.contactPointA,
      rj = this.contactPointB,
      t = this.t,
      G = this.G;
    // G = [-t -rixt t rjxt]
    // And remember, this is a pure velocity constraint, g is always zero!
    G[0] = -t[0];
    G[1] = -t[1];
    G[2] = -vec2.crossLength(ri, t);
    G[3] = t[0];
    G[4] = t[1];
    G[5] = vec2.crossLength(rj, t);
    let GW = this.computeGW(),
      GiMf = this.computeGiMf();
    let B = /* - g * a  */ - GW * b - h * GiMf;
    return B;
  }
}, { "../math/vec2": 30, "../utils/Utils": 57, "./Equation": 22 }], 24: [function (_dereq_, module, exports) {
  let Equation = _dereq_("./Equation"),
    vec2 = _dereq_('../math/vec2');
  module.exports = RotationalLockEquation;
  /**Locks the relative angle between two bodies. The constraint tries to keep the dot product between two vectors, local in each body, to zero. The local angle in body i is a parameter.
  *
  * @class RotationalLockEquation
  * @constructor
  * @extends Equation
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {Object} [options]
  * @param {Number} [options.angle] Angle to add to the local vector in bodyA.
  */
  function RotationalLockEquation(bodyA, bodyB, options) {
    options = options || {}
    Equation.call(this, bodyA, bodyB, -Number.MAX_VALUE, Number.MAX_VALUE);
    /**@property {number} angle
    */
    this.angle = options.angle || 0;
    let G = this.G;
    G[2] = 1;
    G[5] = -1;
  }
  RotationalLockEquation.prototype = new Equation();
  RotationalLockEquation.prototype.constructor = RotationalLockEquation;
  let worldVectorA = vec2.create(),
    worldVectorB = vec2.create(),
    xAxis = vec2.fromValues(1, 0),
    yAxis = vec2.fromValues(0, 1);
  RotationalLockEquation.prototype.computeGq() {
    vec2.rotate(worldVectorA, xAxis, this.bodyA.angle + this.angle);
    vec2.rotate(worldVectorB, yAxis, this.bodyB.angle);
    return vec2.dot(worldVectorA, worldVectorB);
  }
}, { "../math/vec2": 30, "./Equation": 22 }], 25: [function (_dereq_, module, exports) {
  let Equation = _dereq_("./Equation"),
    vec2 = _dereq_('../math/vec2');
  module.exports = RotationalVelocityEquation;
  /**Syncs rotational velocity of two bodies, or sets a relative velocity (motor).
  *
  * @class RotationalVelocityEquation
  * @constructor
  * @extends Equation
  * @param {Body} bodyA
  * @param {Body} bodyB
  */
  function RotationalVelocityEquation(bodyA, bodyB) {
    Equation.call(this, bodyA, bodyB, -Number.MAX_VALUE, Number.MAX_VALUE);
    this.relativeVelocity = 1;
    this.ratio = 1;
  }
  RotationalVelocityEquation.prototype = new Equation();
  RotationalVelocityEquation.prototype.constructor = RotationalVelocityEquation;
  RotationalVelocityEquation.prototype.computeB(a, b, h) {
    let G = this.G;
    G[2] = -1;
    G[5] = this.ratio;
    let GiMf = this.computeGiMf();
    let GW = this.computeGW();
    let B = - GW * b - h * GiMf;
    return B;
  }
}, { "../math/vec2": 30, "./Equation": 22 }], 26: [function (_dereq_, module, exports) {
  /**Base class for objects that dispatches events.
  * @class EventEmitter
  * @constructor
  */
  let EventEmitter () { }
  module.exports = EventEmitter;
  EventEmitter.prototype = {
    constructor: EventEmitter,
    /**Add an event listener
    * @param {String} type
    * @param {Function} listener
    * @returns {EventEmitter} The self object, for chainability.
    */
    on: function (type, listener, context) {
      listener.context = context || this;
      if (this._listeners === undefined) {
        this._listeners = {}
      }
      let listeners = this._listeners;
      if (listeners[type] === undefined) {
        listeners[type] = [];
      }
      if (listeners[type].indexOf(listener) === - 1) {
        listeners[type].push(listener);
      }
      return this;
    },
    /**Check if an event listener is added
    * @param {String} type
    * @param {Function} listener
    * @returns {Boolean}
    */
    has: function (type, listener) {
      if (this._listeners === undefined) {
        return false;
      }
      let listeners = this._listeners;
      if (listener) {
        if (listeners[type] !== undefined && listeners[type].indexOf(listener) !== - 1) {
          return true;
        }
      } else {
        if (listeners[type] !== undefined) {
          return true;
        }
      }
      return false;
    },
    /**Remove an event listener
    * @param {String} type
    * @param {Function} listener
    * @returns {EventEmitter} The self object, for chainability.
    */
    off: function (type, listener) {
      if (this._listeners === undefined) {
        return this;
      }
      let listeners = this._listeners;
      let index = listeners[type].indexOf(listener);
      if (index !== - 1) {
        listeners[type].splice(index, 1);
      }
      return this;
    },
    /**Emit an event.
    * @param {Object} event
    * @param {String} event.type
    * @returns {EventEmitter} The self object, for chainability.
    */
    emit: function (event) {
      if (this._listeners === undefined) {
        return this;
      }
      let listeners = this._listeners;
      let listenerArray = listeners[event.type];
      if (listenerArray !== undefined) {
        event.target = this;
        for (let i = 0, l = listenerArray.length; i < l; i++) {
          let listener = listenerArray[i];
          listener.call(listener.context, event);
        }
      }
      return this;
    }
  }
}, {}], 27: [function (_dereq_, module, exports) {
  let Material = _dereq_('./Material');
  let Equation = _dereq_('../equations/Equation');
  module.exports = ContactMaterial;
  /**Defines what happens when two materials meet, such as what friction coefficient to use. You can also set other things such as restitution, surface velocity and constraint parameters.
  * @class ContactMaterial
  * @constructor
  * @param {Material} materialA
  * @param {Material} materialB
  * @param {Object}   [options]
  * @param {Number}   [options.friction=0.3]       Friction coefficient.
  * @param {Number}   [options.restitution=0]      Restitution coefficient aka "bounciness".
  * @param {Number}   [options.stiffness]          ContactEquation stiffness.
  * @param {Number}   [options.relaxation]         ContactEquation relaxation.
  * @param {Number}   [options.frictionStiffness]  FrictionEquation stiffness.
  * @param {Number}   [options.frictionRelaxation] FrictionEquation relaxation.
  * @param {Number}   [options.surfaceVelocity=0]  Surface velocity.
  * @author schteppe
  */
  function ContactMaterial(materialA, materialB, options) {
    options = options || {}
    if (!(materialA instanceof Material) || !(materialB instanceof Material)) {
      throw new Error("First two arguments must be Material instances.");
    }
    /**The contact material identifier
    * @property id
    * @type {Number}
    */
    this.id = ContactMaterial.idCounter++;
    /**First material participating in the contact material
    * @property materialA
    * @type {Material}
    */
    this.materialA = materialA;
    /**Second material participating in the contact material
    * @property materialB
    * @type {Material}
    */
    this.materialB = materialB;
    /**Friction coefficient to use in the contact of these two materials. Friction = 0 will make the involved objects super slippery, and friction = 1 will make it much less slippery. A friction coefficient larger than 1 will allow for very large friction forces, which can be convenient for preventing car tires not slip on the ground.
    * @property friction
    * @type {Number}
    * @default 0.3
    */
    this.friction = typeof (options.friction) !== "undefined" ? Number(options.friction) : 0.3;
    /**Restitution, or "bounciness" to use in the contact of these two materials. A restitution of 0 will make no bounce, while restitution=1 will approximately bounce back with the same velocity the object came with.
    * @property restitution
    * @type {Number}
    * @default 0
    */
    this.restitution = typeof (options.restitution) !== "undefined" ? Number(options.restitution) : 0;
    /**Hardness of the contact. Less stiffness will make the objects penetrate more, and will make the contact act more like a spring than a contact force. Default value is {{#crossLink "Equation/DEFAULT_STIFFNESS:property"}}Equation.DEFAULT_STIFFNESS{{/crossLink}}.
    * @property stiffness
    * @type {Number}
    */
    this.stiffness = typeof (options.stiffness) !== "undefined" ? Number(options.stiffness) : Equation.DEFAULT_STIFFNESS;
    /**Relaxation of the resulting ContactEquation that this ContactMaterial generate. Default value is {{#crossLink "Equation/DEFAULT_RELAXATION:property"}}Equation.DEFAULT_RELAXATION{{/crossLink}}.
    * @property relaxation
    * @type {Number}
    */
    this.relaxation = typeof (options.relaxation) !== "undefined" ? Number(options.relaxation) : Equation.DEFAULT_RELAXATION;
    /**Stiffness of the resulting friction force. For most cases, the value of this property should be a large number. I cannot think of any case where you would want less frictionStiffness. Default value is {{#crossLink "Equation/DEFAULT_STIFFNESS:property"}}Equation.DEFAULT_STIFFNESS{{/crossLink}}.
    * @property frictionStiffness
    * @type {Number}
    */
    this.frictionStiffness = typeof (options.frictionStiffness) !== "undefined" ? Number(options.frictionStiffness) : Equation.DEFAULT_STIFFNESS;
    /**Relaxation of the resulting friction force. The default value should be good for most simulations. Default value is {{#crossLink "Equation/DEFAULT_RELAXATION:property"}}Equation.DEFAULT_RELAXATION{{/crossLink}}.
    * @property frictionRelaxation
    * @type {Number}
    */
    this.frictionRelaxation = typeof (options.frictionRelaxation) !== "undefined" ? Number(options.frictionRelaxation) : Equation.DEFAULT_RELAXATION;
    /**Will add surface velocity to this material. If bodyA rests on top if bodyB, and the surface velocity is positive, bodyA will slide to the right.
    * @property {Number} surfaceVelocity
    * @default 0
    */
    this.surfaceVelocity = typeof (options.surfaceVelocity) !== "undefined" ? Number(options.surfaceVelocity) : 0;
    /**Offset to be set on ContactEquations. A positive value will make the bodies penetrate more into each other. Can be useful in scenes where contacts need to be more persistent, for example when stacking. Aka "cure for nervous contacts".
    * @property contactSkinSize
    * @type {Number}
    */
    this.contactSkinSize = 0.005;
  }
  ContactMaterial.idCounter = 0;
}, { "../equations/Equation": 22, "./Material": 28 }], 28: [function (_dereq_, module, exports) {
  module.exports = Material;
  /**Defines a physics material.
  * @class Material
  * @constructor
  * @param {number} id Material identifier
  * @author schteppe
  */
  function Material(id) {
    /**The material identifier
    * @property id
    * @type {Number}
    */
    this.id = id || Material.idCounter++;
  }
  Material.idCounter = 0;
}, {}], 29: [function (_dereq_, module, exports) {
  /*
      PolyK library
      url: http://polyk.ivank.net
      Released under MIT licence.
      Copyright (c) 2012 Ivan Kuckir
      Permission is hereby granted, free of charge, to any person
      obtaining a copy of this software and associated documentation
      files (the "Software"), to deal in the Software without
      restriction, including without limitation the rights to use,
      copy, modify, merge, publish, distribute, sublicense, and/or sell
      copies of the Software, and to permit persons to whom the
      Software is furnished to do so, subject to the following
      conditions:
      The above copyright notice and this permission notice shall be
      included in all copies or substantial portions of the Software.
      THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
      EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
      OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
      NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
      HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
      WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
      FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
      OTHER DEALINGS IN THE SOFTWARE.
  */
  let PolyK = {}
  /*
      Is Polygon self-intersecting?
      O(n^2)
  */
  /*
  PolyK.IsSimple = function(p)
  {
      let n = p.length>>1;
      if(n<4) return true;
      let a1 = new PolyK._P(), a2 = new PolyK._P();
      let b1 = new PolyK._P(), b2 = new PolyK._P();
      let c = new PolyK._P();
      for(let i=0; i<n; i++)
      {
          a1.x = p[2*i  ];
          a1.y = p[2*i+1];
          if(i==n-1)  { a2.x = p[0    ];  a2.y = p[1    ]; }
          else        { a2.x = p[2*i+2];  a2.y = p[2*i+3]; }
          for(let j=0; j<n; j++)
          {
              if(Math.abs(i-j) < 2) continue;
              if(j==n-1 && i==0) continue;
              if(i==n-1 && j==0) continue;
              b1.x = p[2*j  ];
              b1.y = p[2*j+1];
              if(j==n-1)  { b2.x = p[0    ];  b2.y = p[1    ]; }
              else        { b2.x = p[2*j+2];  b2.y = p[2*j+3]; }
              if(PolyK._GetLineIntersection(a1,a2,b1,b2,c) != null) return false;
          }
      }
      return true;
  }
  PolyK.IsConvex = function(p)
  {
      if(p.length<6) return true;
      let l = p.length - 4;
      for(let i=0; i<l; i+=2)
          if(!PolyK._convex(p[i], p[i+1], p[i+2], p[i+3], p[i+4], p[i+5])) return false;
      if(!PolyK._convex(p[l  ], p[l+1], p[l+2], p[l+3], p[0], p[1])) return false;
      if(!PolyK._convex(p[l+2], p[l+3], p[0  ], p[1  ], p[2], p[3])) return false;
      return true;
  }
  */
  PolyK.GetArea(p) {
    if (p.length < 6) return 0;
    let l = p.length - 2;
    let sum = 0;
    for (let i = 0; i < l; i += 2)
      sum += (p[i + 2] - p[i]) * (p[i + 1] + p[i + 3]);
    sum += (p[0] - p[l]) * (p[l + 1] + p[1]);
    return - sum * 0.5;
  }
  /*
  PolyK.GetAABB = function(p)
  {
      let minx = Infinity;
      let miny = Infinity;
      let maxx = -minx;
      let maxy = -miny;
      for(let i=0; i<p.length; i+=2)
      {
          minx = Math.min(minx, p[i  ]);
          maxx = Math.max(maxx, p[i  ]);
          miny = Math.min(miny, p[i+1]);
          maxy = Math.max(maxy, p[i+1]);
      }
      return {x:minx, y:miny, width:maxx-minx, height:maxy-miny}
  }
  */
  PolyK.Triangulate(p) {
    let n = p.length >> 1;
    if (n < 3) return [];
    let tgs = [];
    let avl = [];
    for (let i = 0; i < n; i++) avl.push(i);
    let i = 0;
    let al = n;
    while (al > 3) {
      let i0 = avl[(i + 0) % al];
      let i1 = avl[(i + 1) % al];
      let i2 = avl[(i + 2) % al];
      let ax = p[2 * i0], ay = p[2 * i0 + 1];
      let bx = p[2 * i1], by = p[2 * i1 + 1];
      let cx = p[2 * i2], cy = p[2 * i2 + 1];
      let earFound = false;
      if (PolyK._convex(ax, ay, bx, by, cx, cy)) {
        earFound = true;
        for (let j = 0; j < al; j++) {
          let vi = avl[j];
          if (vi == i0 || vi == i1 || vi == i2) continue;
          if (PolyK._PointInTriangle(p[2 * vi], p[2 * vi + 1], ax, ay, bx, by, cx, cy)) { earFound = false; break; }
        }
      }
      if (earFound) {
        tgs.push(i0, i1, i2);
        avl.splice((i + 1) % al, 1);
        al--;
        i = 0;
      }
      else if (i++ > 3 * al) break;      // no convex angles :(
    }
    tgs.push(avl[0], avl[1], avl[2]);
    return tgs;
  }
  /*
  PolyK.ContainsPoint = function(p, px, py)
  {
      let n = p.length>>1;
      let ax, ay, bx = p[2*n-2]-px, by = p[2*n-1]-py;
      let depth = 0;
      for(let i=0; i<n; i++)
      {
          ax = bx;  ay = by;
          bx = p[2*i  ] - px;
          by = p[2*i+1] - py;
          if(ay< 0 && by< 0) continue;    // both "up" or both "donw"
          if(ay>=0 && by>=0) continue;    // both "up" or both "donw"
          if(ax< 0 && bx< 0) continue;
          let lx = ax + (bx-ax)*(-ay)/(by-ay);
          if(lx>0) depth++;
      }
      return (depth & 1) == 1;
  }
  PolyK.Slice = function(p, ax, ay, bx, by)
  {
      if(PolyK.ContainsPoint(p, ax, ay) || PolyK.ContainsPoint(p, bx, by)) return [p.slice(0)];
      let a = new PolyK._P(ax, ay);
      let b = new PolyK._P(bx, by);
      let iscs = [];  // intersections
      let ps = [];    // points
      for(let i=0; i<p.length; i+=2) ps.push(new PolyK._P(p[i], p[i+1]));
      for(let i=0; i<ps.length; i++)
      {
          let isc = new PolyK._P(0,0);
          isc = PolyK._GetLineIntersection(a, b, ps[i], ps[(i+1)%ps.length], isc);
          if(isc)
          {
              isc.flag = true;
              iscs.push(isc);
              ps.splice(i+1,0,isc);
              i++;
          }
      }
      if(iscs.length == 0) return [p.slice(0)];
      let comp = function(u,v) {return PolyK._P.dist(a,u) - PolyK._P.dist(a,v); }
      iscs.sort(comp);
      let pgs = [];
      let dir = 0;
      while(iscs.length > 0)
      {
          let n = ps.length;
          let i0 = iscs[0];
          let i1 = iscs[1];
          let ind0 = ps.indexOf(i0);
          let ind1 = ps.indexOf(i1);
          let solved = false;
          if(PolyK._firstWithFlag(ps, ind0) == ind1) solved = true;
          else
          {
              i0 = iscs[1];
              i1 = iscs[0];
              ind0 = ps.indexOf(i0);
              ind1 = ps.indexOf(i1);
              if(PolyK._firstWithFlag(ps, ind0) == ind1) solved = true;
          }
          if(solved)
          {
              dir--;
              let pgn = PolyK._getPoints(ps, ind0, ind1);
              pgs.push(pgn);
              ps = PolyK._getPoints(ps, ind1, ind0);
              i0.flag = i1.flag = false;
              iscs.splice(0,2);
              if(iscs.length == 0) pgs.push(ps);
          }
          else { dir++; iscs.reverse(); }
          if(dir>1) break;
      }
      let result = [];
      for(let i=0; i<pgs.length; i++)
      {
          let pg = pgs[i];
          let npg = [];
          for(let j=0; j<pg.length; j++) npg.push(pg[j].x, pg[j].y);
          result.push(npg);
      }
      return result;
  }
  PolyK.Raycast = function(p, x, y, dx, dy, isc)
  {
      let l = p.length - 2;
      let tp = PolyK._tp;
      let a1 = tp[0], a2 = tp[1],
      b1 = tp[2], b2 = tp[3], c = tp[4];
      a1.x = x; a1.y = y;
      a2.x = x+dx; a2.y = y+dy;
      if(isc==null) isc = {dist:0, edge:0, norm:{x:0, y:0}, refl:{x:0, y:0}}
      isc.dist = Infinity;
      for(let i=0; i<l; i+=2)
      {
          b1.x = p[i  ];  b1.y = p[i+1];
          b2.x = p[i+2];  b2.y = p[i+3];
          let nisc = PolyK._RayLineIntersection(a1, a2, b1, b2, c);
          if(nisc) PolyK._updateISC(dx, dy, a1, b1, b2, c, i/2, isc);
      }
      b1.x = b2.x;  b1.y = b2.y;
      b2.x = p[0];  b2.y = p[1];
      let nisc = PolyK._RayLineIntersection(a1, a2, b1, b2, c);
      if(nisc) PolyK._updateISC(dx, dy, a1, b1, b2, c, p.length/2, isc);
      return (isc.dist != Infinity) ? isc : null;
  }
  PolyK.ClosestEdge = function(p, x, y, isc)
  {
      let l = p.length - 2;
      let tp = PolyK._tp;
      let a1 = tp[0],
      b1 = tp[2], b2 = tp[3], c = tp[4];
      a1.x = x; a1.y = y;
      if(isc==null) isc = {dist:0, edge:0, point:{x:0, y:0}, norm:{x:0, y:0}}
      isc.dist = Infinity;
      for(let i=0; i<l; i+=2)
      {
          b1.x = p[i  ];  b1.y = p[i+1];
          b2.x = p[i+2];  b2.y = p[i+3];
          PolyK._pointLineDist(a1, b1, b2, i>>1, isc);
      }
      b1.x = b2.x;  b1.y = b2.y;
      b2.x = p[0];  b2.y = p[1];
      PolyK._pointLineDist(a1, b1, b2, l>>1, isc);
      let idst = 1/isc.dist;
      isc.norm.x = (x-isc.point.x)*idst;
      isc.norm.y = (y-isc.point.y)*idst;
      return isc;
  }
  PolyK._pointLineDist = function(p, a, b, edge, isc)
  {
      let x = p.x, y = p.y, x1 = a.x, y1 = a.y, x2 = b.x, y2 = b.y;
      let A = x - x1;
      let B = y - y1;
      let C = x2 - x1;
      let D = y2 - y1;
      let dot = A * C + B * D;
      let len_sq = C * C + D * D;
      let param = dot / len_sq;
      let xx, yy;
      if (param < 0 || (x1 == x2 && y1 == y2)) {
          xx = x1;
          yy = y1;
      }
      else if (param > 1) {
          xx = x2;
          yy = y2;
      }
      else {
          xx = x1 + param * C;
          yy = y1 + param * D;
      }
      let dx = x - xx;
      let dy = y - yy;
      let dst = Math.sqrt(dx * dx + dy * dy);
      if(dst<isc.dist)
      {
          isc.dist = dst;
          isc.edge = edge;
          isc.point.x = xx;
          isc.point.y = yy;
      }
  }
  PolyK._updateISC = function(dx, dy, a1, b1, b2, c, edge, isc)
  {
      let nrl = PolyK._P.dist(a1, c);
      if(nrl<isc.dist)
      {
          let ibl = 1/PolyK._P.dist(b1, b2);
          let nx = -(b2.y-b1.y)*ibl;
          let ny =  (b2.x-b1.x)*ibl;
          let ddot = 2*(dx*nx+dy*ny);
          isc.dist = nrl;
          isc.norm.x = nx;
          isc.norm.y = ny;
          isc.refl.x = -ddot*nx+dx;
          isc.refl.y = -ddot*ny+dy;
          isc.edge = edge;
      }
  }
  PolyK._getPoints = function(ps, ind0, ind1)
  {
      let n = ps.length;
      let nps = [];
      if(ind1<ind0) ind1 += n;
      for(let i=ind0; i<= ind1; i++) nps.push(ps[i%n]);
      return nps;
  }
  PolyK._firstWithFlag = function(ps, ind)
  {
      let n = ps.length;
      while(true)
      {
          ind = (ind+1)%n;
          if(ps[ind].flag) return ind;
      }
  }
  */
  PolyK._PointInTriangle(px, py, ax, ay, bx, by, cx, cy) {
    let Ray.v0x = cx - ax;
  let Ray.v0y = cy - ay;
  let v1x = bx - ax;
  let v1y = by - ay;
  let v2x = px - ax;
  let v2y = py - ay;
  let dot00 = Ray.v0x * Ray.v0x + Ray.v0y * Ray.v0y;
  let dot01 = Ray.v0x * v1x + Ray.v0y * v1y;
  let dot02 = Ray.v0x * v2x + Ray.v0y * v2y;
  let dot11 = v1x * v1x + v1y * v1y;
  let dot12 = v1x * v2x + v1y * v2y;
  let invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
  let u = (dot11 * dot02 - dot01 * dot12) * invDenom;
  let v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    // Check if point is in triangle
    return (u >= 0) && (v >= 0) && (u + v < 1);
  }
/*
PolyK._RayLineIntersection = function(a1, a2, b1, b2, c)
{
    let dax = (a1.x-a2.x), dbx = (b1.x-b2.x);
    let day = (a1.y-a2.y), dby = (b1.y-b2.y);
    let Den = dax*dby - day*dbx;
    if (Den == 0) return null;  // parallel
    let A = (a1.x * a2.y - a1.y * a2.x);
    let B = (b1.x * b2.y - b1.y * b2.x);
    let I = c;
    let iDen = 1/Den;
    I.x = ( A*dbx - dax*B ) * iDen;
    I.y = ( A*dby - day*B ) * iDen;
    if(!PolyK._InRect(I, b1, b2)) return null;
    if((day>0 && I.y>a1.y) || (day<0 && I.y<a1.y)) return null;
    if((dax>0 && I.x>a1.x) || (dax<0 && I.x<a1.x)) return null;
    return I;
}
PolyK._GetLineIntersection = function(a1, a2, b1, b2, c)
{
    let dax = (a1.x-a2.x), dbx = (b1.x-b2.x);
    let day = (a1.y-a2.y), dby = (b1.y-b2.y);
    let Den = dax*dby - day*dbx;
    if (Den == 0) return null;  // parallel
    let A = (a1.x * a2.y - a1.y * a2.x);
    let B = (b1.x * b2.y - b1.y * b2.x);
    let I = c;
    I.x = ( A*dbx - dax*B ) / Den;
    I.y = ( A*dby - day*B ) / Den;
    if(PolyK._InRect(I, a1, a2) && PolyK._InRect(I, b1, b2)) return I;
    return null;
}
PolyK._InRect = function(a, b, c)
{
    if  (b.x == c.x) return (a.y>=Math.min(b.y, c.y) && a.y<=Math.max(b.y, c.y));
    if  (b.y == c.y) return (a.x>=Math.min(b.x, c.x) && a.x<=Math.max(b.x, c.x));
    if(a.x >= Math.min(b.x, c.x) && a.x <= Math.max(b.x, c.x)
    && a.y >= Math.min(b.y, c.y) && a.y <= Math.max(b.y, c.y))
    return true;
    return false;
}
*/
PolyK._convex(ax, ay, bx, by, cx, cy) {
  return (ay - by) * (cx - bx) + (bx - ax) * (cy - by) >= 0;
}
/*
PolyK._P = function(x,y)
{
    this.x = x;
    this.y = y;
    this.flag = false;
}
PolyK._P.prototype.toString = function()
{
    return "Point ["+this.x+", "+this.y+"]";
}
PolyK._P.dist = function(a,b)
{
    let dx = b.x-a.x;
    let dy = b.y-a.y;
    return Math.sqrt(dx*dx + dy*dy);
}
PolyK._tp = [];
for(let i=0; i<10; i++) PolyK._tp.push(new PolyK._P(0,0));
    */
module.exports = PolyK;
}, { }], 30: [function (_dereq_, module, exports) {
  /* Copyright (c) 2013, Brandon Jones, Colin MacKenzie IV. All rights reserved.
  Redistribution and use in source and binary forms, with or without modification,
  are permitted provided that the following conditions are met:
    * Redistributions of source code must retain the above copyright notice, this
      list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.
  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR
  ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
  ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. */
  /**The vec2 object from glMatrix, with some extensions and some removed methods. See http://glmatrix.net.
  * @class vec2
  */
  let vec2 = module.exports = {}
  let Utils = _dereq_('../utils/Utils');
  /**Make a cross product and only return the z component
   * @param {Array} a
  * @param {Array} b
  * @returns {Number}
  */
  vec2.crossLength(a, b) {
    return a[0] * b[1] - a[1] * b[0];
  }
  /**Cross product between a vector and the Z component of a vector
   * @param {Array} out
  * @param {Array} vec
  * @param {Number} zcomp
  * @returns {Number}
  */
  vec2.crossVZ(out, vec, zcomp) {
    vec2.rotate(out, vec, -Math.PI / 2);// Rotate according to the right hand rule
    vec2.scale(out, out, zcomp);      // Scale with z
    return out;
  }
  /**Cross product between a vector and the Z component of a vector
   * @param {Array} out
  * @param {Number} zcomp
  * @param {Array} vec
  * @returns {Number}
  */
  vec2.crossZV(out, zcomp, vec) {
    vec2.rotate(out, vec, Math.PI / 2); // Rotate according to the right hand rule
    vec2.scale(out, out, zcomp);      // Scale with z
    return out;
  }
  /**Rotate a vector by an angle
   * @param {Array} out
  * @param {Array} a
  * @param {Number} angle
  */
  vec2.rotate(out, a, angle) {
    if (angle !== 0) {
      let c = Math.cos(angle),
        s = Math.sin(angle),
        x = a[0],
        y = a[1];
      out[0] = c * x - s * y;
      out[1] = s * x + c * y;
    } else {
      out[0] = a[0];
      out[1] = a[1];
    }
  }
  /**Rotate a vector 90 degrees clockwise
   * @param {Array} out
  * @param {Array} a
  * @param {Number} angle
  */
  vec2.rotate90cw(out, a) {
    let x = a[0];
    let y = a[1];
    out[0] = y;
    out[1] = -x;
  }
  /**Transform a point position to local frame.
  * @param {Array} out
  * @param {Array} worldPoint
  * @param {Array} framePosition
  * @param {Number} frameAngle
  */
  vec2.toLocalFrame(out, worldPoint, framePosition, frameAngle) {
    vec2.copy(out, worldPoint);
    vec2.sub(out, out, framePosition);
    vec2.rotate(out, out, -frameAngle);
  }
  /**Transform a point position to global frame.
  * @param {Array} out
  * @param {Array} localPoint
  * @param {Array} framePosition
  * @param {Number} frameAngle
  */
  vec2.toGlobalFrame(out, localPoint, framePosition, frameAngle) {
    vec2.copy(out, localPoint);
    vec2.rotate(out, out, frameAngle);
    vec2.add(out, out, framePosition);
  }
  /**Transform a vector to local frame.
  * @param {Array} out
  * @param {Array} worldVector
  * @param {Number} frameAngle
  */
  vec2.vectorToLocalFrame(out, worldVector, frameAngle) {
    vec2.rotate(out, worldVector, -frameAngle);
  }
  /**Transform a point position to global frame.
  * @param {Array} out
  * @param {Array} localVector
  * @param {Number} frameAngle
  */
  vec2.vectorToGlobalFrame(out, localVector, frameAngle) {
    vec2.rotate(out, localVector, frameAngle);
  }
  /**Compute centroid of a triangle spanned by vectors a,b,c. See http://easycalculation.com/analytical/learn-centroid.php
   * @param {Array} out
  * @param {Array} a
  * @param {Array} b
  * @param {Array} c
  * @returns {Array} The out object
  */
  vec2.centroid(out, a, b, c) {
    vec2.add(out, a, b);
    vec2.add(out, out, c);
    vec2.scale(out, out, 1 / 3);
    return out;
  }
  /**Creates a new, empty vec2
   * @returns {Array} a new 2D vector
  */
  vec2.create() {
    let out = new Utils.ARRAY_TYPE(2);
    out[0] = 0;
    out[1] = 0;
    return out;
  }
  /**Creates a new vec2 initialized with values from an existing vector
   * @param {Array} a vector to clone
  * @returns {Array} a new 2D vector
  */
  vec2.clone(a) {
    let out = new Utils.ARRAY_TYPE(2);
    out[0] = a[0];
    out[1] = a[1];
    return out;
  }
  /**Creates a new vec2 initialized with the given values
   * @param {Number} x X component
  * @param {Number} y Y component
  * @returns {Array} a new 2D vector
  */
  vec2.fromValues(x, y) {
    let out = new Utils.ARRAY_TYPE(2);
    out[0] = x;
    out[1] = y;
    return out;
  }
  /**Copy the values from one vec2 to another
   * @param {Array} out the receiving vector
  * @param {Array} a the source vector
  * @returns {Array} out
  */
  vec2.copy(out, a) {
    out[0] = a[0];
    out[1] = a[1];
    return out;
  }
  /**Set the components of a vec2 to the given values
   * @param {Array} out the receiving vector
  * @param {Number} x X component
  * @param {Number} y Y component
  * @returns {Array} out
  */
  vec2.set(out, x, y) {
    out[0] = x;
    out[1] = y;
    return out;
  }
  /**Adds two vec2's
   * @param {Array} out the receiving vector
  * @param {Array} a the first operand
  * @param {Array} b the second operand
  * @returns {Array} out
  */
  vec2.add(out, a, b) {
    out[0] = a[0] + b[0];
    out[1] = a[1] + b[1];
    return out;
  }
  /**Subtracts two vec2's
   * @param {Array} out the receiving vector
  * @param {Array} a the first operand
  * @param {Array} b the second operand
  * @returns {Array} out
  */
  vec2.subtract(out, a, b) {
    out[0] = a[0] - b[0];
    out[1] = a[1] - b[1];
    return out;
  }
  /**Alias for vec2.subtract
   */
  vec2.sub = vec2.subtract;
  /**Multiplies two vec2's
   * @param {Array} out the receiving vector
  * @param {Array} a the first operand
  * @param {Array} b the second operand
  * @returns {Array} out
  */
  vec2.multiply(out, a, b) {
    out[0] = a[0] * b[0];
    out[1] = a[1] * b[1];
    return out;
  }
  /**Alias for vec2.multiply
   */
  vec2.mul = vec2.multiply;
  /**Divides two vec2's
   * @param {Array} out the receiving vector
  * @param {Array} a the first operand
  * @param {Array} b the second operand
  * @returns {Array} out
  */
  vec2.divide(out, a, b) {
    out[0] = a[0] / b[0];
    out[1] = a[1] / b[1];
    return out;
  }
  /**Alias for vec2.divide
   */
  vec2.div = vec2.divide;
  /**Scales a vec2 by a scalar number
   * @param {Array} out the receiving vector
  * @param {Array} a the vector to scale
  * @param {Number} b amount to scale the vector by
  * @returns {Array} out
  */
  vec2.scale(out, a, b) {
    out[0] = a[0] * b;
    out[1] = a[1] * b;
    return out;
  }
  /**Calculates the euclidian distance between two vec2's
   * @param {Array} a the first operand
  * @param {Array} b the second operand
  * @returns {Number} distance between a and b
  */
  vec2.distance(a, b) {
    let x = b[0] - a[0],
      y = b[1] - a[1];
    return Math.sqrt(x * x + y * y);
  }
  /**Alias for vec2.distance
   */
  vec2.dist = vec2.distance;
  /**Calculates the squared euclidian distance between two vec2's
   * @param {Array} a the first operand
  * @param {Array} b the second operand
  * @returns {Number} squared distance between a and b
  */
  vec2.squaredDistance(a, b) {
    let x = b[0] - a[0],
      y = b[1] - a[1];
    return x * x + y * y;
  }
  /**Alias for vec2.squaredDistance
   */
  vec2.sqrDist = vec2.squaredDistance;
  /**Calculates the length of a vec2
   * @param {Array} a vector to calculate length of
  * @returns {Number} length of a
  */
  vec2.length(a) {
    let x = a[0],
      y = a[1];
    return Math.sqrt(x * x + y * y);
  }
  /**Alias for vec2.length
   */
  vec2.len = vec2.length;
  /**Calculates the squared length of a vec2
   * @param {Array} a vector to calculate squared length of
  * @returns {Number} squared length of a
  */
  vec2.squaredLength(a) {
    let x = a[0],
      y = a[1];
    return x * x + y * y;
  }
  /**Alias for vec2.squaredLength
   */
  vec2.sqrLen = vec2.squaredLength;
  /**Negates the components of a vec2
   * @param {Array} out the receiving vector
  * @param {Array} a vector to negate
  * @returns {Array} out
  */
  vec2.negate(out, a) {
    out[0] = -a[0];
    out[1] = -a[1];
    return out;
  }
  /**Normalize a vec2
   * @param {Array} out the receiving vector
  * @param {Array} a vector to normalize
  * @returns {Array} out
  */
  vec2.normalize(out, a) {
    let x = a[0],
      y = a[1];
    let len = x * x + y * y;
    if (len > 0) {
      //TODO: evaluate use of glm_invsqrt here?
      len = 1 / Math.sqrt(len);
      out[0] = a[0] * len;
      out[1] = a[1] * len;
    }
    return out;
  }
  /**Calculates the dot product of two vec2's
   * @param {Array} a the first operand
  * @param {Array} b the second operand
  * @returns {Number} dot product of a and b
  */
  vec2.dot(a, b) {
    return a[0] * b[0] + a[1] * b[1];
  }
  /**Returns a string representation of a vector
   * @param {Array} vec vector to represent as a string
  * @returns {String} string representation of the vector
  */
  vec2.str(a) {
    return 'vec2(' + a[0] + ', ' + a[1] + ')';
  }
  /**Linearly interpolate/mix two vectors.
   * @param {Array} out
  * @param {Array} a First vector
  * @param {Array} b Second vector
  * @param {number} t Lerp factor
  */
  vec2.lerp(out, a, b, t) {
    let ax = a[0],
      ay = a[1];
    out[0] = ax + t * (b[0] - ax);
    out[1] = ay + t * (b[1] - ay);
    return out;
  }
  /**Reflect a vector along a normal.
   * @param {Array} out
  * @param {Array} vector
  * @param {Array} normal
  */
  vec2.reflect(out, vector, normal) {
    let dot = vector[0] * normal[0] + vector[1] * normal[1];
    out[0] = vector[0] - 2 * normal[0] * dot;
    out[1] = vector[1] - 2 * normal[1] * dot;
  }
  /**Get the intersection point between two line segments.
   * @param {Array} out
  * @param {Array} p0
  * @param {Array} p1
  * @param {Array} p2
  * @param {Array} p3
  * @returns {Boolean} True if there was an intersection, otherwise false.
  */
  vec2.getLineSegmentsIntersection(out, p0, p1, p2, p3) {
    let t = vec2.getLineSegmentsIntersectionFraction(p0, p1, p2, p3);
    if (t < 0) {
      return false;
    } else {
      out[0] = p0[0] + (t * (p1[0] - p0[0]));
      out[1] = p0[1] + (t * (p1[1] - p0[1]));
      return true;
    }
  }
  /**(p1 - p0)
   * @param {Array} p0
  * @param {Array} p1
  * @param {Array} p2
  * @param {Array} p3
  * @returns {number} A number between 0 and 1 if there was an intersection, otherwise -1.
  */
  vec2.getLineSegmentsIntersectionFraction(p0, p1, p2, p3) {
    let s1_x = p1[0] - p0[0];
    let s1_y = p1[1] - p0[1];
    let s2_x = p3[0] - p2[0];
    let s2_y = p3[1] - p2[1];
    let s, t;
    s = (-s1_y * (p0[0] - p2[0]) + s1_x * (p0[1] - p2[1])) / (-s2_x * s1_y + s1_x * s2_y);
    t = (s2_x * (p0[1] - p2[1]) - s2_y * (p0[0] - p2[0])) / (-s2_x * s1_y + s1_x * s2_y);
    if (s >= 0 && s <= 1 && t >= 0 && t <= 1) { // Collision detected
      return t;
    }
    return -1; // No collision
  }
}, { "../utils/Utils": 57 }], 31: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2')
    , decomp = _dereq_('poly-decomp')
    , Convex = _dereq_('../shapes/Convex')
    , RaycastResult = _dereq_('../collision/RaycastResult')
    , Ray = _dereq_('../collision/Ray')
    , AABB = _dereq_('../collision/AABB')
    , EventEmitter = _dereq_('../events/EventEmitter');
  module.exports = Body;
  /**A rigid body. Has got a center of mass, position, velocity and a number of
  * shapes that are used for collisions.
  *
  * @class Body
  * @constructor
  * @extends EventEmitter
  * @param {Array} [options.force]
  * @param {Array} [options.position]
  * @param {Array} [options.velocity]
  * @param {Boolean} [options.allowSleep]
  * @param {Boolean} [options.collisionResponse]
  * @param {Number} [options.angle=0]
  * @param {Number} [options.angularForce=0]
  * @param {Number} [options.angularVelocity=0]
  * @param {Number} [options.ccdIterations=10]
  * @param {Number} [options.ccdSpeedThreshold=-1]
  * @param {Number} [options.fixedRotation=false]
  * @param {Number} [options.gravityScale]
  * @param {Number} [options.id]
  * @param {Number} [options.mass=0] A number >= 0. If zero, the .type will be set to Body.STATIC.
  * @param {Number} [options.sleepSpeedLimit]
  * @param {Number} [options.sleepTimeLimit]
  * @param {Object} [options]
  *
  * @example
  *
  *     // Create a typical dynamic body
  *     let body = new Body({
  *         mass: 1,
  *         position: [0, 0],
  *         angle: 0,
  *         velocity: [0, 0],
  *         angularVelocity: 0
  *     });
  *
  *     // Add a circular shape to the body
  *     body.addShape(new Circle({ radius: 1 }));
  *
  *     // Add the body to the world
  *     world.addBody(body);
  */
  function Body(options) {
    options = options || {}
    EventEmitter.call(this);
    /**The body identifyer
    * @property id
    * @type {Number}
    */
    this.id = options.id || ++Body._idCounter;
    /**The world that this body is added to. This property is set to NULL if the body is not added to any world.
    * @property world
    * @type {World}
    */
    this.world = null;
    /**The shapes of the body.
    *
    * @property shapes
    * @type {Array}
    */
    this.shapes = [];
    /**The mass of the body.
    * @property mass
    * @type {number}
    */
    this.mass = options.mass || 0;
    /**The inverse mass of the body.
    * @property invMass
    * @type {number}
    */
    this.invMass = 0;
    /**The inertia of the body around the Z axis.
    * @property inertia
    * @type {number}
    */
    this.inertia = 0;
    /**The inverse inertia of the body.
    * @property invInertia
    * @type {number}
    */
    this.invInertia = 0;
    this.invMassSolve = 0;
    this.invInertiaSolve = 0;
    /**Set to true if you want to fix the rotation of the body.
    * @property fixedRotation
    * @type {Boolean}
    */
    this.fixedRotation = !!options.fixedRotation;
    /**Set to true if you want to fix the body movement along the X axis. The body will still be able to move along Y.
    * @property {Boolean} fixedX
    */
    this.fixedX = !!options.fixedX;
    /**Set to true if you want to fix the body movement along the Y axis. The body will still be able to move along X.
    * @property {Boolean} fixedY
    */
    this.fixedY = !!options.fixedY;
    /**@private
    * @property {array} massMultiplier
    */
    this.massMultiplier = vec2.create();
    /**The position of the body
    * @property position
    * @type {Array}
    */
    this.position = vec2.fromValues(0, 0);
    if (options.position) {
      vec2.copy(this.position, options.position);
    }
    /**The interpolated position of the body. Use this for rendering.
    * @property interpolatedPosition
    * @type {Array}
    */
    this.interpolatedPosition = vec2.fromValues(0, 0);
    /**The interpolated angle of the body. Use this for rendering.
    * @property interpolatedAngle
    * @type {Number}
    */
    this.interpolatedAngle = 0;
    /**The previous position of the body.
    * @property previousPosition
    * @type {Array}
    */
    this.previousPosition = vec2.fromValues(0, 0);
    /**The previous angle of the body.
    * @property previousAngle
    * @type {Number}
    */
    this.previousAngle = 0;
    /**The current velocity of the body.
    * @property velocity
    * @type {Array}
    */
    this.velocity = vec2.fromValues(0, 0);
    if (options.velocity) {
      vec2.copy(this.velocity, options.velocity);
    }
    /**Constraint velocity that was added to the body during the last step.
    * @property vlambda
    * @type {Array}
    */
    this.vlambda = vec2.fromValues(0, 0);
    /**Angular constraint velocity that was added to the body during last step.
    * @property wlambda
    * @type {Array}
    */
    this.wlambda = 0;
    /**The angle of the body, in radians.
    * @property angle
    * @type {number}
    * @example
    *     // The angle property is not normalized to the interval 0 to 2*pi, it can be any value.
    *     // If you need a value between 0 and 2*pi, use the following function to normalize it.
    *     function normalizeAngle(angle){
    *         angle = angle % (2*Math.PI);
    *         if(angle < 0){
    *             angle += (2*Math.PI);
    *         }
    *         return angle;
    *     }
    */
    this.angle = options.angle || 0;
    /**The angular velocity of the body, in radians per second.
    * @property angularVelocity
    * @type {number}
    */
    this.angularVelocity = options.angularVelocity || 0;
    /**The force acting on the body. Since the body force (and {{#crossLink "Body/angularForce:property"}}{{/crossLink}}) will be zeroed after each step, so you need to set the force before each step.
    * @property force
    * @type {Array}
    *
    * @example
    *     // This produces a forcefield of 1 Newton in the positive x direction.
    *     for(let i=0; i<numSteps; i++){
    *         body.force[0] = 1;
    *         world.step(1/60);
    *     }
    *
    * @example
    *     // This will apply a rotational force on the body
    *     for(let i=0; i<numSteps; i++){
    *         body.angularForce = -3;
    *         world.step(1/60);
    *     }
    */
    this.force = vec2.create();
    if (options.force) {
      vec2.copy(this.force, options.force);
    }
    /**The angular force acting on the body. See {{#crossLink "Body/force:property"}}{{/crossLink}}.
    * @property angularForce
    * @type {number}
    */
    this.angularForce = options.angularForce || 0;
    /**The linear damping acting on the body in the velocity direction. Should be a value between 0 and 1.
    * @property damping
    * @type {Number}
    * @default 0.1
    */
    this.damping = typeof (options.damping) === "number" ? options.damping : 0.1;
    /**The angular force acting on the body. Should be a value between 0 and 1.
    * @property angularDamping
    * @type {Number}
    * @default 0.1
    */
    this.angularDamping = typeof (options.angularDamping) === "number" ? options.angularDamping : 0.1;
    /**The type of motion this body has. Should be one of: {{#crossLink "Body/STATIC:property"}}Body.STATIC{{/crossLink}}, {{#crossLink "Body/DYNAMIC:property"}}Body.DYNAMIC{{/crossLink}} and {{#crossLink "Body/KINEMATIC:property"}}Body.KINEMATIC{{/crossLink}}.
    *
    * * Static bodies do not move, and they do not respond to forces or collision.
    * * Dynamic bodies body can move and respond to collisions and forces.
    * * Kinematic bodies only moves according to its .velocity, and does not respond to collisions or force.
    *
    * @property type
    * @type {number}
    *
    * @example
    *     // Bodies are static by default. Static bodies will never move.
    *     let body = new Body();
    *     console.log(body.type == Body.STATIC); // true
    *
    * @example
    *     // By setting the mass of a body to a nonzero number, the body
    *     // will become dynamic and will move and interact with other bodies.
    *     let dynamicBody = new Body({
    *         mass : 1
    *     });
    *     console.log(dynamicBody.type == Body.DYNAMIC); // true
    *
    * @example
    *     // Kinematic bodies will only move if you change their velocity.
    *     let kinematicBody = new Body({
    *         type: Body.KINEMATIC // Type can be set via the options object.
    *     });
    */
    this.type = Body.STATIC;
    if (typeof (options.type) !== 'undefined') {
      this.type = options.type;
    } else if (!options.mass) {
      this.type = Body.STATIC;
    } else {
      this.type = Body.DYNAMIC;
    }
    /**Bounding circle radius.
    * @property boundingRadius
    * @type {Number}
    */
    this.boundingRadius = 0;
    /**Bounding box of this body.
    * @property aabb
    * @type {AABB}
    */
    this.aabb = new AABB();
    /**Indicates if the AABB needs update. Update it with {{#crossLink "Body/updateAABB:method"}}.updateAABB(){{/crossLink}}.
    * @property aabbNeedsUpdate
    * @type {Boolean}
    * @see updateAABB
    *
    * @example
    *     // Force update the AABB
    *     body.aabbNeedsUpdate = true;
    *     body.updateAABB();
    *     console.log(body.aabbNeedsUpdate); // false
    */
    this.aabbNeedsUpdate = true;
    /**If true, the body will automatically fall to sleep. Note that you need to enable sleeping in the {{#crossLink "World"}}{{/crossLink}} before anything will happen.
    * @property allowSleep
    * @type {Boolean}
    * @default true
    */
    this.allowSleep = options.allowSleep !== undefined ? options.allowSleep : true;
    this.wantsToSleep = false;
    /**One of {{#crossLink "Body/AWAKE:property"}}Body.AWAKE{{/crossLink}}, {{#crossLink "Body/SLEEPY:property"}}Body.SLEEPY{{/crossLink}} and {{#crossLink "Body/SLEEPING:property"}}Body.SLEEPING{{/crossLink}}.
    *
    * The body is initially Body.AWAKE. If its velocity norm is below .sleepSpeedLimit, the sleepState will become Body.SLEEPY. If the body continues to be Body.SLEEPY for .sleepTimeLimit seconds, it will fall asleep (Body.SLEEPY).
    *
    * @property sleepState
    * @type {Number}
    * @default Body.AWAKE
    */
    this.sleepState = Body.AWAKE;
    /**If the speed (the norm of the velocity) is smaller than this value, the body is considered sleepy.
    * @property sleepSpeedLimit
    * @type {Number}
    * @default 0.2
    */
    this.sleepSpeedLimit = options.sleepSpeedLimit !== undefined ? options.sleepSpeedLimit : 0.2;
    /**If the body has been sleepy for this sleepTimeLimit seconds, it is considered sleeping.
    * @property sleepTimeLimit
    * @type {Number}
    * @default 1
    */
    this.sleepTimeLimit = options.sleepTimeLimit !== undefined ? options.sleepTimeLimit : 1;
    /**Gravity scaling factor. If you want the body to ignore gravity, set this to zero. If you want to reverse gravity, set it to -1.
    * @property {Number} gravityScale
    * @default 1
    */
    this.gravityScale = options.gravityScale !== undefined ? options.gravityScale : 1;
    /**Whether to produce contact forces when in contact with other bodies. Note that contacts will be generated, but they will be disabled. That means that this body will move through other bodies, but it will still trigger contact events, etc.
    * @property {Boolean} collisionResponse
    */
    this.collisionResponse = options.collisionResponse !== undefined ? options.collisionResponse : true;
    /**How long the body has been sleeping.
    * @property {Number} idleTime
    */
    this.idleTime = 0;
    /**The last time when the body went to SLEEPY state.
    * @property {Number} timeLastSleepy
    * @private
    */
    this.timeLastSleepy = 0;
    /**If the body speed exceeds this threshold, CCD (continuous collision detection) will be enabled. Set it to a negative number to disable CCD completely for this body.
    * @property {number} ccdSpeedThreshold
    * @default -1
    */
    this.ccdSpeedThreshold = options.ccdSpeedThreshold !== undefined ? options.ccdSpeedThreshold : -1;
    /**The number of iterations that should be used when searching for the time of impact during CCD. A larger number will assure that there's a small penetration on CCD collision, but a small number will give more performance.
    * @property {number} ccdIterations
    * @default 10
    */
    this.ccdIterations = options.ccdIterations !== undefined ? options.ccdIterations : 10;
    this.concavePath = null;
    this._wakeUpAfterNarrowphase = false;
    this.updateMassProperties();
  }
  Body.prototype = new EventEmitter();
  Body.prototype.constructor = Body;
  Body._idCounter = 0;
  /**@private
  */
  Body.prototype.updateSolveMassProperties() {
    if (this.sleepState === Body.SLEEPING || this.type === Body.KINEMATIC) {
      this.invMassSolve = 0;
      this.invInertiaSolve = 0;
    } else {
      this.invMassSolve = this.invMass;
      this.invInertiaSolve = this.invInertia;
    }
  }
  /**Set the total density of the body
  * @param {number} density
  */
  Body.prototype.setDensity(density) {
    let totalArea = this.getArea();
    this.mass = totalArea * density;
    this.updateMassProperties();
  }
  /**Get the total area of all shapes in the body
  * @returns {Number}
  */
  Body.prototype.getArea() {
    let totalArea = 0;
    for (let i = 0; i < this.shapes.length; i++) {
      totalArea += this.shapes[i].area;
    }
    return totalArea;
  }
  /**Get the AABB from the body. The AABB is updated if necessary.
  * @returns {AABB} The AABB instance (this.aabb)
  */
  Body.prototype.getAABB() {
    if (this.aabbNeedsUpdate) {
      this.updateAABB();
    }
    return this.aabb;
  }
  let shapeAABB = new AABB(),
    tmp = vec2.create();
  /**Updates the AABB of the Body, and set .aabbNeedsUpdate = false.
  */
  Body.prototype.updateAABB() {
    let shapes = this.shapes,
      N = shapes.length,
      offset = tmp,
      bodyAngle = this.angle;
    for (let i = 0; i !== N; i++) {
      let shape = shapes[i],
        angle = shape.angle + bodyAngle;
      // Get shape world offset
      vec2.rotate(offset, shape.position, bodyAngle);
      vec2.add(offset, offset, this.position);
      // Get shape AABB
      shape.computeAABB(shapeAABB, offset, angle);
      if (i === 0) {
        this.aabb.copy(shapeAABB);
      } else {
        this.aabb.extend(shapeAABB);
      }
    }
    this.aabbNeedsUpdate = false;
  }
  /**Update the bounding radius of the body (this.boundingRadius). Should be done if any of the shape dimensions or positions are changed.
  */
  Body.prototype.updateBoundingRadius() {
    let shapes = this.shapes,
      N = shapes.length,
      radius = 0;
    for (let i = 0; i !== N; i++) {
      let shape = shapes[i],
        offset = vec2.length(shape.position),
        r = shape.boundingRadius;
      if (offset + r > radius) {
        radius = offset + r;
      }
    }
    this.boundingRadius = radius;
  }
  /**Add a shape to the body. You can pass a local transform when adding a shape,
  * so that the shape gets an offset and angle relative to the body center of mass.
  * Will automatically update the mass properties and bounding radius.
  *
  * @param {Shape}              shape
  * @param {Array} [offset] Local body offset of the shape.
  * @param {Number}             [angle]  Local body angle.
  *
  * @example
  *     let body = new Body(),
  *         shape = new Circle({ radius: 1 });
  *
  *     // Add the shape to the body, positioned in the center
  *     body.addShape(shape);
  *
  *     // Add another shape to the body, positioned 1 unit length from the body center of mass along the local x-axis.
  *     body.addShape(shape,[1,0]);
  *
  *     // Add another shape to the body, positioned 1 unit length from the body center of mass along the local y-axis, and rotated 90 degrees CCW.
  *     body.addShape(shape,[0,1],Math.PI/2);
  */
  Body.prototype.addShape(shape, offset, angle) {
    if (shape.body) {
      throw new Error('A shape can only be added to one body.');
    }
    shape.body = this;
    // Copy the offset vector
    if (offset) {
      vec2.copy(shape.position, offset);
    } else {
      vec2.set(shape.position, 0, 0);
    }
    shape.angle = angle || 0;
    this.shapes.push(shape);
    this.updateMassProperties();
    this.updateBoundingRadius();
    this.aabbNeedsUpdate = true;
  }
  /**Remove a shape
  * @param {Shape} shape
  * @returns {Boolean} True if the shape was found and removed, else false.
  */
  Body.prototype.removeShape(shape) {
    let idx = this.shapes.indexOf(shape);
    if (idx !== -1) {
      this.shapes.splice(idx, 1);
      this.aabbNeedsUpdate = true;
      shape.body = null;
      return true;
    } else {
      return false;
    }
  }
  /**Updates .inertia, .invMass, .invInertia for this Body. Should be called when
  * changing the structure or mass of the Body.
  *
  *
  * @example
  *     body.mass += 1;
  *     body.updateMassProperties();
  */
  Body.prototype.updateMassProperties() {
    if (this.type === Body.STATIC || this.type === Body.KINEMATIC) {
      this.mass = Number.MAX_VALUE;
      this.invMass = 0;
      this.inertia = Number.MAX_VALUE;
      this.invInertia = 0;
    } else {
      let shapes = this.shapes,
        N = shapes.length,
        m = this.mass / N,
        I = 0;
      if (!this.fixedRotation) {
        for (let i = 0; i < N; i++) {
          let shape = shapes[i],
            r2 = vec2.squaredLength(shape.position),
            Icm = shape.computeMomentOfInertia(m);
          I += Icm + m * r2;
        }
        this.inertia = I;
        this.invInertia = I > 0 ? 1 / I : 0;
      } else {
        this.inertia = Number.MAX_VALUE;
        this.invInertia = 0;
      }
      // Inverse mass properties are easy
      this.invMass = 1 / this.mass;
      vec2.set(
        this.massMultiplier,
        this.fixedX ? 0 : 1,
        this.fixedY ? 0 : 1
      );
    }
  }
  let Body_applyForce_r = vec2.create();
  /**Apply force to a point relative to the center of mass of the body. This could for example be a point on the RigidBody surface. Applying force this way will add to Body.force and Body.angularForce. If relativePoint is zero, the force will be applied directly on the center of mass, and the torque produced will be zero.
  * @param {Array} force The force to add.
  * @param {Array} [relativePoint] A world point to apply the force on.
  */
  Body.prototype.applyForce(force, relativePoint) {
    // Add linear force
    vec2.add(this.force, this.force, force);
    if (relativePoint) {
      // Compute produced rotational force
      let rotForce = vec2.crossLength(relativePoint, force);
      // Add rotational force
      this.angularForce += rotForce;
    }
  }
  /**Apply force to a body-local point.
  * @param {Array} localForce The force vector to add, oriented in local body space.
  * @param {Array} [localPoint] A point relative to the body in world space. If not given, it is set to zero and all of the impulse will be excerted on the center of mass.
  */
  let Body_applyForce_forceWorld = vec2.create();
  let Body_applyForce_pointWorld = vec2.create();
  let Body_applyForce_pointLocal = vec2.create();
  Body.prototype.applyForceLocal(localForce, localPoint) {
    localPoint = localPoint || Body_applyForce_pointLocal;
    let worldForce = Body_applyForce_forceWorld;
    let worldPoint = Body_applyForce_pointWorld;
    this.vectorToWorldFrame(worldForce, localForce);
    this.vectorToWorldFrame(worldPoint, localPoint);
    this.applyForce(worldForce, worldPoint);
  }
  /**time). Impulses will be added to Body.velocity and Body.angularVelocity.
  * @param {Array} impulse The impulse vector to add, oriented in world space.
  * @param {Array} [relativePoint] A point relative to the body in world space. If not given, it is set to zero and all of the impulse will be excerted on the center of mass.
  */
  let Body_applyImpulse_velo = vec2.create();
  Body.prototype.applyImpulse(impulseVector, relativePoint) {
    if (this.type !== Body.DYNAMIC) {
      return;
    }
    // Compute produced central impulse velocity
    let velo = Body_applyImpulse_velo;
    vec2.scale(velo, impulseVector, this.invMass);
    vec2.multiply(velo, this.massMultiplier, velo);
    // Add linear impulse
    vec2.add(this.velocity, velo, this.velocity);
    if (relativePoint) {
      // Compute produced rotational impulse velocity
      let rotVelo = vec2.crossLength(relativePoint, impulseVector);
      rotVelo *= this.invInertia;
      // Add rotational Impulse
      this.angularVelocity += rotVelo;
    }
  }
  /**time). Impulses will be added to Body.velocity and Body.angularVelocity.
  * @param {Array} impulse The impulse vector to add, oriented in world space.
  * @param {Array} [relativePoint] A point relative to the body in world space. If not given, it is set to zero and all of the impulse will be excerted on the center of mass.
  */
  let Body_applyImpulse_impulseWorld = vec2.create();
  let Body_applyImpulse_pointWorld = vec2.create();
  let Body_applyImpulse_pointLocal = vec2.create();
  Body.prototype.applyImpulseLocal(localImpulse, localPoint) {
    localPoint = localPoint || Body_applyImpulse_pointLocal;
    let worldImpulse = Body_applyImpulse_impulseWorld;
    let worldPoint = Body_applyImpulse_pointWorld;
    this.vectorToWorldFrame(worldImpulse, localImpulse);
    this.vectorToWorldFrame(worldPoint, localPoint);
    this.applyImpulse(worldImpulse, worldPoint);
  }
  /**Transform a world point to local body frame.
  * @param {Array} out          The vector to store the result in
  * @param {Array} worldPoint   The input world point
  */
  Body.prototype.toLocalFrame(out, worldPoint) {
    vec2.toLocalFrame(out, worldPoint, this.position, this.angle);
  }
  /**Transform a local point to world frame.
  * @param {Array} out          The vector to store the result in
  * @param {Array} localPoint   The input local point
  */
  Body.prototype.toWorldFrame(out, localPoint) {
    vec2.toGlobalFrame(out, localPoint, this.position, this.angle);
  }
  /**Transform a world point to local body frame.
  * @param {Array} out          The vector to store the result in
  * @param {Array} worldVector  The input world vector
  */
  Body.prototype.vectorToLocalFrame(out, worldVector) {
    vec2.vectorToLocalFrame(out, worldVector, this.angle);
  }
  /**Transform a local point to world frame.
  * @param {Array} out          The vector to store the result in
  * @param {Array} localVector  The input local vector
  */
  Body.prototype.vectorToWorldFrame(out, localVector) {
    vec2.vectorToGlobalFrame(out, localVector, this.angle);
  }
  /**Reads a polygon shape path, and assembles convex shapes from that and puts them at proper offset points.
  * @param {Array} path An array of 2d vectors, e.g. [[0,0],[0,1],...] that resembles a concave or convex polygon. The shape must be simple and without holes.
  * @param {Object} [options]
  * @param {Boolean} [options.optimalDecomp=false]   Set to true if you need optimal decomposition. Warning: very slow for polygons with more than 10 vertices.
  * @param {Boolean} [options.skipSimpleCheck=false] Set to true if you already know that the path is not intersecting itself.
  * @param {Boolean|Number} [options.removeCollinearPoints=false] Set to a number (angle threshold value) to remove collinear points, or false to keep all points.
  * @returns {Boolean} True on success, else false.
  */
  Body.prototype.fromPolygon(path, options) {
    options = options || {}
    // Remove all shapes
    for (let i = this.shapes.length; i >= 0; --i) {
      this.removeShape(this.shapes[i]);
    }
    let p = new decomp.Polygon();
    p.vertices = path;
    // Make it counter-clockwise
    p.makeCCW();
    if (typeof (options.removeCollinearPoints) === "number") {
      p.removeCollinearPoints(options.removeCollinearPoints);
    }
    // Check if any line segment intersects the path itself
    if (typeof (options.skipSimpleCheck) === "undefined") {
      if (!p.isSimple()) {
        return false;
      }
    }
    // Save this path for later
    this.concavePath = p.vertices.slice(0);
    for (let i = 0; i < this.concavePath.length; i++) {
      let v = [0, 0];
      vec2.copy(v, this.concavePath[i]);
      this.concavePath[i] = v;
    }
    // Slow or fast decomp?
    let convexes;
    if (options.optimalDecomp) {
      convexes = p.decomp();
    } else {
      convexes = p.quickDecomp();
    }
    let cm = vec2.create();
    // Add convexes
    for (let i = 0; i !== convexes.length; i++) {
      // Create convex
      let c = new Convex({ vertices: convexes[i].vertices });
      // Move all vertices so its center of mass is in the local center of the convex
      for (let j = 0; j !== c.vertices.length; j++) {
        let v = c.vertices[j];
        vec2.sub(v, v, c.centerOfMass);
      }
      vec2.scale(cm, c.centerOfMass, 1);
      c.updateTriangles();
      c.updateCenterOfMass();
      c.updateBoundingRadius();
      // Add the shape
      this.addShape(c, cm);
    }
    this.adjustCenterOfMass();
    this.aabbNeedsUpdate = true;
    return true;
  }
  let adjustCenterOfMass_tmp1 = vec2.fromValues(0, 0),
    adjustCenterOfMass_tmp2 = vec2.fromValues(0, 0),
    adjustCenterOfMass_tmp3 = vec2.fromValues(0, 0),
    adjustCenterOfMass_tmp4 = vec2.fromValues(0, 0);
  /**Moves the shape offsets so their center of mass becomes the body center of mass.
  */
  Body.prototype.adjustCenterOfMass() {
    let offset_times_area = adjustCenterOfMass_tmp2,
      sum = adjustCenterOfMass_tmp3,
      cm = adjustCenterOfMass_tmp4,
      totalArea = 0;
    vec2.set(sum, 0, 0);
    for (let i = 0; i !== this.shapes.length; i++) {
      let s = this.shapes[i];
      vec2.scale(offset_times_area, s.position, s.area);
      vec2.add(sum, sum, offset_times_area);
      totalArea += s.area;
    }
    vec2.scale(cm, sum, 1 / totalArea);
    // Now move all shapes
    for (let i = 0; i !== this.shapes.length; i++) {
      let s = this.shapes[i];
      vec2.sub(s.position, s.position, cm);
    }
    // Move the body position too
    vec2.add(this.position, this.position, cm);
    // And concave path
    for (let i = 0; this.concavePath && i < this.concavePath.length; i++) {
      vec2.sub(this.concavePath[i], this.concavePath[i], cm);
    }
    this.updateMassProperties();
    this.updateBoundingRadius();
  }
  /**Sets the force on the body to zero.
  */
  Body.prototype.setZeroForce() {
    vec2.set(this.force, 0.0, 0.0);
    this.angularForce = 0.0;
  }
  Body.prototype.resetConstraintVelocity() {
    let b = this,
      vlambda = b.vlambda;
    vec2.set(vlambda, 0, 0);
    b.wlambda = 0;
  }
  Body.prototype.addConstraintVelocity() {
    let b = this,
      v = b.velocity;
    vec2.add(v, v, b.vlambda);
    b.angularVelocity += b.wlambda;
  }
  /**Apply damping, see <a href="http://code.google.com/p/bullet/issues/detail?id=74">this</a> for details.
  * @param {number} dt Current time step
  */
  Body.prototype.applyDamping(dt) {
    if (this.type === Body.DYNAMIC) { // Only for dynamic bodies
      let v = this.velocity;
      vec2.scale(v, v, Math.pow(1.0 - this.damping, dt));
      this.angularVelocity *= Math.pow(1.0 - this.angularDamping, dt);
    }
  }
  /**Wake the body up. Normally you should not need this, as the body is automatically awoken at events such as collisions.
  * Sets the sleepState to {{#crossLink "Body/AWAKE:property"}}Body.AWAKE{{/crossLink}} and emits the wakeUp event if the body wasn't awake before.
  */
  Body.prototype.wakeUp() {
    let s = this.sleepState;
    this.sleepState = Body.AWAKE;
    this.idleTime = 0;
    if (s !== Body.AWAKE) {
      this.emit(Body.wakeUpEvent);
    }
  }
  /**Force body sleep
  */
  Body.prototype.sleep() {
    this.sleepState = Body.SLEEPING;
    this.angularVelocity = 0;
    this.angularForce = 0;
    vec2.set(this.velocity, 0, 0);
    vec2.set(this.force, 0, 0);
    this.emit(Body.sleepEvent);
  }
  /**Called every timestep to update internal sleep timer and change sleep state if needed.
  * @param {number} time The world time in seconds
  * @param {Boolean} dontSleep
  * @param {number} dt
  */
  Body.prototype.sleepTick(time, dontSleep, dt) {
    if (!this.allowSleep || this.type === Body.SLEEPING) {
      return;
    }
    this.wantsToSleep = false;
    let sleepState = this.sleepState,
      speedSquared = vec2.squaredLength(this.velocity) + Math.pow(this.angularVelocity, 2),
      speedLimitSquared = Math.pow(this.sleepSpeedLimit, 2);
    // Add to idle time
    if (speedSquared >= speedLimitSquared) {
      this.idleTime = 0;
      this.sleepState = Body.AWAKE;
    } else {
      this.idleTime += dt;
      this.sleepState = Body.SLEEPY;
    }
    if (this.idleTime > this.sleepTimeLimit) {
      if (!dontSleep) {
        this.sleep();
      } else {
        this.wantsToSleep = true;
      }
    }
  }
  /**Check if the body is overlapping another body. Note that this method only works if the body was added to a World and if at least one step was taken.
  * @param {Body} body
  * @returns {Boolean}
  */
  Body.prototype.overlaps(body) {
    return this.world.overlapKeeper.bodiesAreOverlapping(this, body);
  }
  let integrate_fhMinv = vec2.create();
  let integrate_velodt = vec2.create();
  /**Move the body forward in time given its current velocity.
  * @param {Number} dt
  */
  Body.prototype.integrate(dt) {
    let minv = this.invMass,
      f = this.force,
      pos = this.position,
      velo = this.velocity;
    // Save old position
    vec2.copy(this.previousPosition, this.position);
    this.previousAngle = this.angle;
    // Velocity update
    if (!this.fixedRotation) {
      this.angularVelocity += this.angularForce * this.invInertia * dt;
    }
    vec2.scale(integrate_fhMinv, f, dt * minv);
    vec2.multiply(integrate_fhMinv, this.massMultiplier, integrate_fhMinv);
    vec2.add(velo, integrate_fhMinv, velo);
    // CCD
    if (!this.integrateToTimeOfImpact(dt)) {
      // Regular position update
      vec2.scale(integrate_velodt, velo, dt);
      vec2.add(pos, pos, integrate_velodt);
      if (!this.fixedRotation) {
        this.angle += this.angularVelocity * dt;
      }
    }
    this.aabbNeedsUpdate = true;
  }
  let result = new RaycastResult();
  let ray = new Ray({
    mode: Ray.ALL
  });
  let direction = vec2.create();
  let end = vec2.create();
  let startToEnd = vec2.create();
  let rememberPosition = vec2.create();
  Body.prototype.integrateToTimeOfImpact(dt) {
    if (this.ccdSpeedThreshold < 0 || vec2.squaredLength(this.velocity) < Math.pow(this.ccdSpeedThreshold, 2)) {
      return false;
    }
    vec2.normalize(direction, this.velocity);
    vec2.scale(end, this.velocity, dt);
    vec2.add(end, end, this.position);
    vec2.sub(startToEnd, end, this.position);
    let startToEndAngle = this.angularVelocity * dt;
    let len = vec2.length(startToEnd);
    let timeOfImpact = 1;
    let hit;
    let that = this;
    result.reset();
    ray.callback(result) {
      if (result.body === that) {
        return;
      }
      hit = result.body;
      result.getHitPoint(end, ray);
      vec2.sub(startToEnd, end, that.position);
      timeOfImpact = vec2.length(startToEnd) / len;
      result.stop();
    }
    vec2.copy(ray.from, this.position);
    vec2.copy(ray.to, end);
    ray.update();
    this.world.raycast(result, ray);
    if (!hit) {
      return false;
    }
    let rememberAngle = this.angle;
    vec2.copy(rememberPosition, this.position);
    // Got a start and end point. Approximate time of impact using binary search
    let iter = 0;
    let tmin = 0;
    let tmid = 0;
    let tmax = timeOfImpact;
    while (tmax >= tmin && iter < this.ccdIterations) {
      iter++;
      // calculate the midpoint
      tmid = (tmax - tmin) / 2;
      // Move the body to that point
      vec2.scale(integrate_velodt, startToEnd, timeOfImpact);
      vec2.add(this.position, rememberPosition, integrate_velodt);
      this.angle = rememberAngle + startToEndAngle * timeOfImpact;
      this.updateAABB();
      // check overlap
      let overlaps = this.aabb.overlaps(hit.aabb) && this.world.narrowphase.bodiesOverlap(this, hit);
      if (overlaps) {
        // change min to search upper interval
        tmin = tmid;
      } else {
        // change max to search lower interval
        tmax = tmid;
      }
    }
    timeOfImpact = tmid;
    vec2.copy(this.position, rememberPosition);
    this.angle = rememberAngle;
    // move to TOI
    vec2.scale(integrate_velodt, startToEnd, timeOfImpact);
    vec2.add(this.position, this.position, integrate_velodt);
    if (!this.fixedRotation) {
      this.angle += startToEndAngle * timeOfImpact;
    }
    return true;
  }
  /**Get velocity of a point in the body.
  * @param {Array} result A vector to store the result in
  * @param {Array} relativePoint A world oriented vector, indicating the position of the point to get the velocity from
  * @returns {Array} The result vector
  */
  Body.prototype.getVelocityAtPoint(result, relativePoint) {
    vec2.crossVZ(result, relativePoint, this.angularVelocity);
    vec2.subtract(result, this.velocity, result);
    return result;
  }
  /**@event sleepy
  */
  Body.sleepyEvent = {
    type: "sleepy"
  }
  /**@event sleep
  */
  Body.sleepEvent = {
    type: "sleep"
  }
  /**@event wakeup
  */
  Body.wakeUpEvent = {
    type: "wakeup"
  }
  /**Dynamic body.
  * @property DYNAMIC
  * @type {Number}
   */
  Body.DYNAMIC = 1;
  /**Static body.
  * @property STATIC
  * @type {Number}
   */
  Body.STATIC = 2;
  /**Kinematic body.
  * @property KINEMATIC
  * @type {Number}
   */
  Body.KINEMATIC = 4;
  /**@property AWAKE
  * @type {Number}
   */
  Body.AWAKE = 0;
  /**@property SLEEPY
  * @type {Number}
   */
  Body.SLEEPY = 1;
  /**@property SLEEPING
  * @type {Number}
   */
  Body.SLEEPING = 2;

}, { "../collision/AABB": 7, "../collision/Ray": 11, "../collision/RaycastResult": 12, "../events/EventEmitter": 26, "../math/vec2": 30, "../shapes/Convex": 40, "poly-decomp": 5 }], 32: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2');
  let Spring = _dereq_('./Spring');
  let Utils = _dereq_('../utils/Utils');
  module.exports = LinearSpring;
  /**A spring, connecting two bodies.
  *
  * The Spring explicitly adds force and angularForce to the bodies.
  *
  * @class LinearSpring
  * @extends Spring
  * @constructor
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {Object} [options]
  * @param {number} [options.restLength]   A number > 0. Default is the current distance between the world anchor points.
  * @param {number} [options.stiffness=100]  Spring constant (see Hookes Law). A number >= 0.
  * @param {number} [options.damping=1]      A number >= 0. Default: 1
  * @param {Array}  [options.worldAnchorA]   Where to hook the spring to body A, in world coordinates. Overrides the option "localAnchorA" if given.
  * @param {Array}  [options.worldAnchorB]
  * @param {Array}  [options.localAnchorA]   Where to hook the spring to body A, in local body coordinates. Defaults to the body center.
  * @param {Array}  [options.localAnchorB]
  */
  function LinearSpring(bodyA, bodyB, options) {
    options = options || {}
    Spring.call(this, bodyA, bodyB, options);
    /**Anchor for bodyA in local bodyA coordinates.
    * @property localAnchorA
    * @type {Array}
    */
    this.localAnchorA = vec2.fromValues(0, 0);
    /**Anchor for bodyB in local bodyB coordinates.
    * @property localAnchorB
    * @type {Array}
    */
    this.localAnchorB = vec2.fromValues(0, 0);
    if (options.localAnchorA) { vec2.copy(this.localAnchorA, options.localAnchorA); }
    if (options.localAnchorB) { vec2.copy(this.localAnchorB, options.localAnchorB); }
    if (options.worldAnchorA) { this.setWorldAnchorA(options.worldAnchorA); }
    if (options.worldAnchorB) { this.setWorldAnchorB(options.worldAnchorB); }
    let worldAnchorA = vec2.create();
    let worldAnchorB = vec2.create();
    this.getWorldAnchorA(worldAnchorA);
    this.getWorldAnchorB(worldAnchorB);
    let worldDistance = vec2.distance(worldAnchorA, worldAnchorB);
    /**Rest length of the spring.
    * @property restLength
    * @type {number}
    */
    this.restLength = typeof (options.restLength) === "number" ? options.restLength : worldDistance;
  }
  LinearSpring.prototype = new Spring();
  LinearSpring.prototype.constructor = LinearSpring;
  /**Set the anchor point on body A, using world coordinates.
  * @param {Array} worldAnchorA
  */
  LinearSpring.prototype.setWorldAnchorA(worldAnchorA) {
    this.bodyA.toLocalFrame(this.localAnchorA, worldAnchorA);
  }
  /**Set the anchor point on body B, using world coordinates.
  * @param {Array} worldAnchorB
  */
  LinearSpring.prototype.setWorldAnchorB(worldAnchorB) {
    this.bodyB.toLocalFrame(this.localAnchorB, worldAnchorB);
  }
  /**Get the anchor point on body A, in world coordinates.
  * @param {Array} result The vector to store the result in.
  */
  LinearSpring.prototype.getWorldAnchorA(result) {
    this.bodyA.toWorldFrame(result, this.localAnchorA);
  }
  /**Get the anchor point on body B, in world coordinates.
  * @param {Array} result The vector to store the result in.
  */
  LinearSpring.prototype.getWorldAnchorB(result) {
    this.bodyB.toWorldFrame(result, this.localAnchorB);
  }
  let applyForce_r = vec2.create(),
    applyForce_r_unit = vec2.create(),
    applyForce_u = vec2.create(),
    applyForce_f = vec2.create(),
    applyForce_worldAnchorA = vec2.create(),
    applyForce_worldAnchorB = vec2.create(),
    applyForce_ri = vec2.create(),
    applyForce_rj = vec2.create(),
    applyForce_tmp = vec2.create();
  /**Apply the spring force to the connected bodies.
  */
  LinearSpring.prototype.applyForce() {
    let k = this.stiffness,
      d = this.damping,
      l = this.restLength,
      bodyA = this.bodyA,
      bodyB = this.bodyB,
      r = applyForce_r,
      r_unit = applyForce_r_unit,
      u = applyForce_u,
      f = applyForce_f,
      tmp = applyForce_tmp;
    let worldAnchorA = applyForce_worldAnchorA,
      worldAnchorB = applyForce_worldAnchorB,
      ri = applyForce_ri,
      rj = applyForce_rj;
    // Get world anchors
    this.getWorldAnchorA(worldAnchorA);
    this.getWorldAnchorB(worldAnchorB);
    // Get offset points
    vec2.sub(ri, worldAnchorA, bodyA.position);
    vec2.sub(rj, worldAnchorB, bodyB.position);
    // Compute distance vector between world anchor points
    vec2.sub(r, worldAnchorB, worldAnchorA);
    let rlen = vec2.len(r);
    vec2.normalize(r_unit, r);
    //console.log(rlen)
    //console.log("A",vec2.str(worldAnchorA),"B",vec2.str(worldAnchorB))
    // Compute relative velocity of the anchor points, u
    vec2.sub(u, bodyB.velocity, bodyA.velocity);
    vec2.crossZV(tmp, bodyB.angularVelocity, rj);
    vec2.add(u, u, tmp);
    vec2.crossZV(tmp, bodyA.angularVelocity, ri);
    vec2.sub(u, u, tmp);
    // F = - k * ( x - L ) - D * ( u )
    vec2.scale(f, r_unit, -k * (rlen - l) - d * vec2.dot(u, r_unit));
    // Add forces to bodies
    vec2.sub(bodyA.force, bodyA.force, f);
    vec2.add(bodyB.force, bodyB.force, f);
    // Angular force
    let ri_x_f = vec2.crossLength(ri, f);
    let rj_x_f = vec2.crossLength(rj, f);
    bodyA.angularForce -= ri_x_f;
    bodyB.angularForce += rj_x_f;
  }
}, { "../math/vec2": 30, "../utils/Utils": 57, "./Spring": 34 }], 33: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2');
  let Spring = _dereq_('./Spring');
  module.exports = RotationalSpring;
  /**A rotational spring, connecting two bodies rotation. This spring explicitly adds angularForce (torque) to the bodies.
  *
  * The spring can be combined with a {{#crossLink "RevoluteConstraint"}}{{/crossLink}} to make, for example, a mouse trap.
  *
  * @class RotationalSpring
  * @extends Spring
  * @constructor
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {Object} [options]
  * @param {number} [options.restAngle] The relative angle of bodies at which the spring is at rest. If not given, it's set to the current relative angle between the bodies.
  * @param {number} [options.stiffness=100] Spring constant (see Hookes Law). A number >= 0.
  * @param {number} [options.damping=1] A number >= 0.
  */
  function RotationalSpring(bodyA, bodyB, options) {
    options = options || {}
    Spring.call(this, bodyA, bodyB, options);
    /**Rest angle of the spring.
    * @property restAngle
    * @type {number}
    */
    this.restAngle = typeof (options.restAngle) === "number" ? options.restAngle : bodyB.angle - bodyA.angle;
  }
  RotationalSpring.prototype = new Spring();
  RotationalSpring.prototype.constructor = RotationalSpring;
  /**Apply the spring force to the connected bodies.
  */
  RotationalSpring.prototype.applyForce() {
    let k = this.stiffness,
      d = this.damping,
      l = this.restAngle,
      bodyA = this.bodyA,
      bodyB = this.bodyB,
      x = bodyB.angle - bodyA.angle,
      u = bodyB.angularVelocity - bodyA.angularVelocity;
    let torque = - k * (x - l) - d * u * 0;
    bodyA.angularForce -= torque;
    bodyB.angularForce += torque;
  }
}, { "../math/vec2": 30, "./Spring": 34 }], 34: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2');
  let Utils = _dereq_('../utils/Utils');
  module.exports = Spring;
  /**A spring, connecting two bodies. The Spring explicitly adds force and angularForce to the bodies and does therefore not put load on the constraint solver.
  *
  * @class Spring
  * @constructor
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @param {Object} [options]
  * @param {number} [options.stiffness=100]  Spring constant (see Hookes Law). A number >= 0.
  * @param {number} [options.damping=1]      A number >= 0. Default: 1
  * @param {Array}  [options.localAnchorA]   Where to hook the spring to body A, in local body coordinates. Defaults to the body center.
  * @param {Array}  [options.localAnchorB]
  * @param {Array}  [options.worldAnchorA]   Where to hook the spring to body A, in world coordinates. Overrides the option "localAnchorA" if given.
  * @param {Array}  [options.worldAnchorB]
  */
  function Spring(bodyA, bodyB, options) {
    options = Utils.defaults(options, {
      stiffness: 100,
      damping: 1,
    });
    /**Stiffness of the spring.
    * @property stiffness
    * @type {number}
    */
    this.stiffness = options.stiffness;
    /**Damping of the spring.
    * @property damping
    * @type {number}
    */
    this.damping = options.damping;
    /**First connected body.
    * @property bodyA
    * @type {Body}
    */
    this.bodyA = bodyA;
    /**Second connected body.
    * @property bodyB
    * @type {Body}
    */
    this.bodyB = bodyB;
  }
  /**Apply the spring force to the connected bodies.
  */
  Spring.prototype.applyForce() {
    // To be implemented by subclasses
  }
}, { "../math/vec2": 30, "../utils/Utils": 57 }], 35: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2');
  let Utils = _dereq_('../utils/Utils');
  let Constraint = _dereq_('../constraints/Constraint');
  let FrictionEquation = _dereq_('../equations/FrictionEquation');
  let Body = _dereq_('../objects/Body');
  module.exports = TopDownVehicle;
  /**@class TopDownVehicle
  * @constructor
  * @param {Body} chassisBody A dynamic body, already added to the world.
  * @param {Object} [options]
  *
  * @example
  *
  *     // Create a dynamic body for the chassis
  *     let chassisBody = new Body({
  *         mass: 1
  *     });
  *     let boxShape = new Box({ width: 0.5, height: 1 });
  *     chassisBody.addShape(boxShape);
  *     world.addBody(chassisBody);
  *
  *     // Create the vehicle
  *     let vehicle = new TopDownVehicle(chassisBody);
  *
  *     // Add one front wheel and one back wheel - we don't actually need four :)
  *     let frontWheel = vehicle.addWheel({
  *         localPosition: [0, 0.5] // front
  *     });
  *     frontWheel.setSideFriction(4);
  *
  *     // Back wheel
  *     let backWheel = vehicle.addWheel({
  *         localPosition: [0, -0.5] // back
  *     });
  *     backWheel.setSideFriction(3); // Less side friction on back wheel makes it easier to drift
  *     vehicle.addToWorld(world);
  *
  *     // Steer value zero means straight forward. Positive is left and negative right.
  *     frontWheel.steerValue = Math.PI / 16;
  *
  *     // Engine force forward
  *     backWheel.engineForce = 10;
  *     backWheel.setBrakeForce(0);
  */
  function TopDownVehicle(chassisBody, options) {
    options = options || {}
    /**@property {Body} chassisBody
    */
    this.chassisBody = chassisBody;
    /**@property {Array} wheels
    */
    this.wheels = [];
    // A dummy body to constrain the chassis to
    this.groundBody = new Body({ mass: 0 });
    this.world = null;
    let that = this;
    this.preStepCallback() {
      that.update();
    }
  }
  /**@param {World} world
  */
  TopDownVehicle.prototype.addToWorld(world) {
    this.world = world;
    world.addBody(this.groundBody);
    world.on('preStep', this.preStepCallback);
    for (let i = 0; i < this.wheels.length; i++) {
      let wheel = this.wheels[i];
      world.addConstraint(wheel);
    }
  }
  /**@param {World} world
  */
  TopDownVehicle.prototype.removeFromWorld() {
    let world = this.world;
    world.removeBody(this.groundBody);
    world.off('preStep', this.preStepCallback);
    for (let i = 0; i < this.wheels.length; i++) {
      let wheel = this.wheels[i];
      world.removeConstraint(wheel);
    }
    this.world = null;
  }
  /**@param {object} [wheelOptions]
  * @returns {WheelConstraint}
  */
  TopDownVehicle.prototype.addWheel(wheelOptions) {
    let wheel = new WheelConstraint(this, wheelOptions);
    this.wheels.push(wheel);
    return wheel;
  }
  TopDownVehicle.prototype.update() {
    for (let i = 0; i < this.wheels.length; i++) {
      this.wheels[i].update();
    }
  }
  /**@class WheelConstraint
  * @constructor
  * @extends {Constraint}
  * @param {Vehicle} vehicle
  * @param {object} [options]
  * @param {Array} [options.localForwardVector]The local wheel forward vector in local body space. Default is zero.
  * @param {Array} [options.localPosition] The local position of the wheen in the chassis body. Default is zero - the center of the body.
  * @param {Array} [options.sideFriction=5] The max friction force in the sideways direction.
  */
  function WheelConstraint(vehicle, options) {
    options = options || {}
    this.vehicle = vehicle;
    this.forwardEquation = new FrictionEquation(vehicle.chassisBody, vehicle.groundBody);
    this.sideEquation = new FrictionEquation(vehicle.chassisBody, vehicle.groundBody);
    /**@property {number} steerValue
    */
    this.steerValue = 0;
    /**@property {number} engineForce
    */
    this.engineForce = 0;
    this.setSideFriction(options.sideFriction !== undefined ? options.sideFriction : 5);
    /**@property {Array} localForwardVector
    */
    this.localForwardVector = vec2.fromValues(0, 1);
    if (options.localForwardVector) {
      vec2.copy(this.localForwardVector, options.localForwardVector);
    }
    /**@property {Array} localPosition
    */
    this.localPosition = vec2.fromValues(0, 0);
    if (options.localPosition) {
      vec2.copy(this.localPosition, options.localPosition);
    }
    Constraint.apply(this, vehicle.chassisBody, vehicle.groundBody);
    this.equations.push(
      this.forwardEquation,
      this.sideEquation
    );
    this.setBrakeForce(0);
  }
  WheelConstraint.prototype = new Constraint();
  WheelsetBrakeForce(force) {
    this.forwardEquation.setSlipForce(force);
  }
  WheelsetSideFriction(force) {
    this.sideEquation.setSlipForce(force);
  }
  let worldVelocity = vec2.create();
  let relativePoint = vec2.create();
  WheelgetSpeed() {
    this.vehicle.chassisBody.vectorToWorldFrame(relativePoint, this.localForwardVector);
    this.vehicle.chassisBody.getVelocityAtPoint(worldVelocity, relativePoint);
    return vec2.dot(worldVelocity, relativePoint);
  }
  let tmpVec = vec2.create();
  Wheelupdate() {
    // Directional
    this.vehicle.chassisBody.vectorToWorldFrame(this.forwardEquation.t, this.localForwardVector);
    vec2.rotate(this.sideEquation.t, this.localForwardVector, Math.PI / 2);
    this.vehicle.chassisBody.vectorToWorldFrame(this.sideEquation.t, this.sideEquation.t);
    vec2.rotate(this.forwardEquation.t, this.forwardEquation.t, this.steerValue);
    vec2.rotate(this.sideEquation.t, this.sideEquation.t, this.steerValue);
    // Attachment point
    this.vehicle.chassisBody.toWorldFrame(this.forwardEquation.contactPointB, this.localPosition);
    vec2.copy(this.sideEquation.contactPointB, this.forwardEquation.contactPointB);
    this.vehicle.chassisBody.vectorToWorldFrame(this.forwardEquation.contactPointA, this.localPosition);
    vec2.copy(this.sideEquation.contactPointA, this.forwardEquation.contactPointA);
    // Add engine force
    vec2.normalize(tmpVec, this.forwardEquation.t);
    vec2.scale(tmpVec, tmpVec, this.engineForce);
    this.vehicle.chassisBody.applyForce(tmpVec, this.forwardEquation.contactPointA);
  }
}, { "../constraints/Constraint": 14, "../equations/FrictionEquation": 23, "../math/vec2": 30, "../objects/Body": 31, "../utils/Utils": 57 }], 36: [function (_dereq_, module, exports) {
  // Export p2 classes
  let p2 = module.exports = {
    AABB: _dereq_('./collision/AABB'),
    AngleLockEquation: _dereq_('./equations/AngleLockEquation'),
    Body: _dereq_('./objects/Body'),
    Broadphase: _dereq_('./collision/Broadphase'),
    Capsule: _dereq_('./shapes/Capsule'),
    Circle: _dereq_('./shapes/Circle'),
    Constraint: _dereq_('./constraints/Constraint'),
    ContactEquation: _dereq_('./equations/ContactEquation'),
    ContactEquationPool: _dereq_('./utils/ContactEquationPool'),
    ContactMaterial: _dereq_('./material/ContactMaterial'),
    Convex: _dereq_('./shapes/Convex'),
    DistanceConstraint: _dereq_('./constraints/DistanceConstraint'),
    Equation: _dereq_('./equations/Equation'),
    EventEmitter: _dereq_('./events/EventEmitter'),
    FrictionEquation: _dereq_('./equations/FrictionEquation'),
    FrictionEquationPool: _dereq_('./utils/FrictionEquationPool'),
    GearConstraint: _dereq_('./constraints/GearConstraint'),
    GSSolver: _dereq_('./solver/GSSolver'),
    Heightfield: _dereq_('./shapes/Heightfield'),
    Line: _dereq_('./shapes/Line'),
    LockConstraint: _dereq_('./constraints/LockConstraint'),
    Material: _dereq_('./material/Material'),
    Narrowphase: _dereq_('./collision/Narrowphase'),
    NaiveBroadphase: _dereq_('./collision/NaiveBroadphase'),
    Particle: _dereq_('./shapes/Particle'),
    Plane: _dereq_('./shapes/Plane'),
    Pool: _dereq_('./utils/Pool'),
    RevoluteConstraint: _dereq_('./constraints/RevoluteConstraint'),
    PrismaticConstraint: _dereq_('./constraints/PrismaticConstraint'),
    Ray: _dereq_('./collision/Ray'),
    RaycastResult: _dereq_('./collision/RaycastResult'),
    Box: _dereq_('./shapes/Box'),
    RotationalVelocityEquation: _dereq_('./equations/RotationalVelocityEquation'),
    SAPBroadphase: _dereq_('./collision/SAPBroadphase'),
    Shape: _dereq_('./shapes/Shape'),
    Solver: _dereq_('./solver/Solver'),
    Spring: _dereq_('./objects/Spring'),
    TopDownVehicle: _dereq_('./objects/TopDownVehicle'),
    LinearSpring: _dereq_('./objects/LinearSpring'),
    RotationalSpring: _dereq_('./objects/RotationalSpring'),
    Utils: _dereq_('./utils/Utils'),
    World: _dereq_('./world/World'),
    vec2: _dereq_('./math/vec2'),
    version: _dereq_('../package.json').version,
  }
  Object.defineProperty(p2, 'Rectangle', {
    get: function () {
      console.warn('The Rectangle class has been renamed to Box.');
      return this.Box;
    }
  });
}, { "../package.json": 6, "./collision/AABB": 7, "./collision/Broadphase": 8, "./collision/NaiveBroadphase": 9, "./collision/Narrowphase": 10, "./collision/Ray": 11, "./collision/RaycastResult": 12, "./collision/SAPBroadphase": 13, "./constraints/Constraint": 14, "./constraints/DistanceConstraint": 15, "./constraints/GearConstraint": 16, "./constraints/LockConstraint": 17, "./constraints/PrismaticConstraint": 18, "./constraints/RevoluteConstraint": 19, "./equations/AngleLockEquation": 20, "./equations/ContactEquation": 21, "./equations/Equation": 22, "./equations/FrictionEquation": 23, "./equations/RotationalVelocityEquation": 25, "./events/EventEmitter": 26, "./material/ContactMaterial": 27, "./material/Material": 28, "./math/vec2": 30, "./objects/Body": 31, "./objects/LinearSpring": 32, "./objects/RotationalSpring": 33, "./objects/Spring": 34, "./objects/TopDownVehicle": 35, "./shapes/Box": 37, "./shapes/Capsule": 38, "./shapes/Circle": 39, "./shapes/Convex": 40, "./shapes/Heightfield": 41, "./shapes/Line": 42, "./shapes/Particle": 43, "./shapes/Plane": 44, "./shapes/Shape": 45, "./solver/GSSolver": 46, "./solver/Solver": 47, "./utils/ContactEquationPool": 48, "./utils/FrictionEquationPool": 49, "./utils/Pool": 55, "./utils/Utils": 57, "./world/World": 61 }], 37: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2')
    , Shape = _dereq_('./Shape')
    , Convex = _dereq_('./Convex');
  module.exports = Box;
  /**Box shape class.
  * @class Box
  * @constructor
  * @param {object} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  * @param {Number} [options.width=1] Total width of the box
  * @param {Number} [options.height=1] Total height of the box
  * @extends Convex
  */
  function Box(options) {
    if (typeof (arguments[0]) === 'number' && typeof (arguments[1]) === 'number') {
      options = {
        width: arguments[0],
        height: arguments[1]
      }
      console.warn('The Rectangle has been renamed to Box and its constructor signature has changed. Please use the following format: new Box({ width: 1, height: 1, ... })');
    }
    options = options || {}
    /**Total width of the box
    * @property width
    * @type {Number}
    */
    let width = this.width = options.width || 1;
    /**Total height of the box
    * @property height
    * @type {Number}
    */
    let height = this.height = options.height || 1;
    let verts = [
      vec2.fromValues(-width / 2, -height / 2),
      vec2.fromValues(width / 2, -height / 2),
      vec2.fromValues(width / 2, height / 2),
      vec2.fromValues(-width / 2, height / 2)
    ];
    let axes = [
      vec2.fromValues(1, 0),
      vec2.fromValues(0, 1)
    ];
    options.vertices = verts;
    options.axes = axes;
    options.type = Shape.BOX;
    Convex.call(this, options);
  }
  Box.prototype = new Convex();
  Box.prototype.constructor = Box;
  /**Compute moment of inertia
  * @param {Number} mass
  * @returns {Number}
  */
  Box.prototype.computeMomentOfInertia(mass) {
    let w = this.width,
      h = this.height;
    return mass * (h * h + w * w) / 12;
  }
  /**Update the bounding radius
  */
  Box.prototype.updateBoundingRadius() {
    let w = this.width,
      h = this.height;
    this.boundingRadius = Math.sqrt(w * w + h * h) / 2;
  }
  let corner1 = vec2.create(),
    corner2 = vec2.create(),
    corner3 = vec2.create(),
    corner4 = vec2.create();
  /**@param {AABB}   out      The resulting AABB.
  * @param {Array}  position
  * @param {Number} angle
  */
  Box.prototype.computeAABB(out, position, angle) {
    out.setFromPoints(this.vertices, position, angle, 0);
  }
  Box.prototype.updateArea() {
    this.area = this.width * this.height;
  }

}, { "../math/vec2": 30, "./Convex": 40, "./Shape": 45 }], 38: [function (_dereq_, module, exports) {
  let Shape = _dereq_('./Shape')
    , vec2 = _dereq_('../math/vec2');
  module.exports = Capsule;
  /**Capsule shape class.
  * @class Capsule
  * @constructor
  * @extends Shape
  * @param {object} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  * @param {Number} [options.length=1] The distance between the end points
  * @param {Number} [options.radius=1] Radius of the capsule
  * @example
  *     let capsuleShape = new Capsule({
  *         length: 1,
  *         radius: 2
  *     });
  *     body.addShape(capsuleShape);
  */
  function Capsule(options) {
    if (typeof (arguments[0]) === 'number' && typeof (arguments[1]) === 'number') {
      options = {
        length: arguments[0],
        radius: arguments[1]
      }
      console.warn('The Capsule constructor signature has changed. Please use the following format: new Capsule({ radius: 1, length: 1 })');
    }
    options = options || {}
    /**The distance between the end points.
    * @property {Number} length
    */
    this.length = options.length || 1;
    /**The radius of the capsule.
    * @property {Number} radius
    */
    this.radius = options.radius || 1;
    options.type = Shape.CAPSULE;
    Shape.call(this, options);
  }
  Capsule.prototype = new Shape();
  Capsule.prototype.constructor = Capsule;
  /**Compute the mass moment of inertia of the Capsule.
  * @param {Number} mass
  * @returns {Number}
  * @todo
  */
  Capsule.prototype.computeMomentOfInertia(mass) {
    // Approximate with rectangle
    let r = this.radius,
      w = this.length + r, // 2*r is too much, 0 is too little
      h = r * 2;
    return mass * (h * h + w * w) / 12;
  }
  Capsule.prototype.updateBoundingRadius() {
    this.boundingRadius = this.radius + this.length / 2;
  }
  Capsule.prototype.updateArea() {
    this.area = Math.PI * this.radius * this.radius + this.radius * 2 * this.length;
  }
  let r = vec2.create();
  /**@param {AABB}   out      The resulting AABB.
  * @param {Array}  position
  * @param {Number} angle
  */
  Capsule.prototype.computeAABB(out, position, angle) {
    let radius = this.radius;
    // Compute center position of one of the the circles, world oriented, but with local offset
    vec2.set(r, this.length / 2, 0);
    if (angle !== 0) {
      vec2.rotate(r, r, angle);
    }
    // Get bounds
    vec2.set(out.upperBound, Math.max(r[0] + radius, -r[0] + radius),
      Math.max(r[1] + radius, -r[1] + radius));
    vec2.set(out.lowerBound, Math.min(r[0] - radius, -r[0] - radius),
      Math.min(r[1] - radius, -r[1] - radius));
    // Add offset
    vec2.add(out.lowerBound, out.lowerBound, position);
    vec2.add(out.upperBound, out.upperBound, position);
  }
  let intersectCapsule_Ray.hitPointWorld = vec2.create();
  let intersectCapsule_normal = vec2.create();
  let intersectCapsule_l0 = vec2.create();
  let intersectCapsule_l1 = vec2.create();
  let intersectCapsule_unit_y = vec2.fromValues(0, 1);
/**@param {RaycastResult} result
* @param {Ray} ray
* @param {array} position
* @param {number} angle
*/
Capsule.prototype.raycast(result, ray, position, angle) {
  let from = ray.from;
  let to = ray.to;
  let direction = ray.direction;
  let Ray.hitPointWorld = intersectCapsule_Ray.hitPointWorld;
  let normal = intersectCapsule_normal;
  let l0 = intersectCapsule_l0;
  let l1 = intersectCapsule_l1;
  // The sides
  let halfLen = this.length / 2;
  for(let i = 0; i < 2; i++) {
    // get start and end of the line
    let y = this.radius * (i * 2 - 1);
    vec2.set(l0, -halfLen, y);
    vec2.set(l1, halfLen, y);
    vec2.toGlobalFrame(l0, l0, position, angle);
    vec2.toGlobalFrame(l1, l1, position, angle);
    let delta = vec2.getLineSegmentsIntersectionFraction(from, to, l0, l1);
    if(delta >= 0) {
  vec2.rotate(normal, intersectCapsule_unit_y, angle);
  vec2.scale(normal, normal, (i * 2 - 1));
  ray.reportIntersection(result, delta, normal, -1);
  if(result.shouldStop(ray)) {
  return;
}
      }
    }
  // Circles
  let diagonalLengthSquared = Math.pow(this.radius, 2) + Math.pow(halfLen, 2);
    for (let i = 0; i < 2; i++) {
  vec2.set(l0, halfLen * (i * 2 - 1), 0);
  vec2.toGlobalFrame(l0, l0, position, angle);
  let a = Math.pow(to[0] - from[0], 2) + Math.pow(to[1] - from[1], 2);
  let b = 2 * ((to[0] - from[0]) * (from[0] - l0[0]) + (to[1] - from[1]) * (from[1] - l0[1]));
  let c = Math.pow(from[0] - l0[0], 2) + Math.pow(from[1] - l0[1], 2) - Math.pow(this.radius, 2);
  let delta = Math.pow(b, 2) - 4 * a * c;
  if (delta < 0) {
    // No intersection
    continue;
  } else if (delta === 0) {
    // single intersection point
    vec2.lerp(Ray.hitPointWorld, from, to, delta);
    if (vec2.squaredDistance(Ray.hitPointWorld, position) > diagonalLengthSquared) {
      vec2.sub(normal, Ray.hitPointWorld, l0);
      vec2.normalize(normal, normal);
      ray.reportIntersection(result, delta, normal, -1);
      if (result.shouldStop(ray)) {
        return;
      }
    }
  } else {
    let sqrtDelta = Math.sqrt(delta);
    let inv2a = 1 / (2 * a);
    let d1 = (- b - sqrtDelta) * inv2a;
    let d2 = (- b + sqrtDelta) * inv2a;
    if (d1 >= 0 && d1 <= 1) {
      vec2.lerp(Ray.hitPointWorld, from, to, d1);
      if (vec2.squaredDistance(Ray.hitPointWorld, position) > diagonalLengthSquared) {
        vec2.sub(normal, Ray.hitPointWorld, l0);
        vec2.normalize(normal, normal);
        ray.reportIntersection(result, d1, normal, -1);
        if (result.shouldStop(ray)) {
          return;
        }
      }
    }
    if (d2 >= 0 && d2 <= 1) {
      vec2.lerp(Ray.hitPointWorld, from, to, d2);
      if (vec2.squaredDistance(Ray.hitPointWorld, position) > diagonalLengthSquared) {
        vec2.sub(normal, Ray.hitPointWorld, l0);
        vec2.normalize(normal, normal);
        ray.reportIntersection(result, d2, normal, -1);
        if (result.shouldStop(ray)) {
          return;
        }
      }
    }
  }
}
  }
}, { "../math/vec2": 30, "./Shape": 45 }], 39: [function (_dereq_, module, exports) {
  let Shape = _dereq_('./Shape')
    , vec2 = _dereq_('../math/vec2');
  module.exports = Circle;
  /**Circle shape class.
  * @class Circle
  * @extends Shape
  * @constructor
  * @param {options} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  * @param {number} [options.radius=1] The radius of this circle
  *
  * @example
  *     let circleShape = new Circle({ radius: 1 });
  *     body.addShape(circleShape);
  */
  function Circle(options) {
    if (typeof (arguments[0]) === 'number') {
      options = {
        radius: arguments[0]
      }
      console.warn('The Circle constructor signature has changed. Please use the following format: new Circle({ radius: 1 })');
    }
    options = options || {}
    /**The radius of the circle.
    * @property radius
    * @type {number}
    */
    this.radius = options.radius || 1;
    options.type = Shape.CIRCLE;
    Shape.call(this, options);
  }
  Circle.prototype = new Shape();
  Circle.prototype.constructor = Circle;
  /**@param {Number} mass
  * @returns {Number}
  */
  Circle.prototype.computeMomentOfInertia(mass) {
    let r = this.radius;
    return mass * r * r / 2;
  }
  /**@returns {Number}
  */
  Circle.prototype.updateBoundingRadius() {
    this.boundingRadius = this.radius;
  }
  /**@returns {Number}
  */
  Circle.prototype.updateArea() {
    this.area = Math.PI * this.radius * this.radius;
  }
  /**@param {AABB}   out      The resulting AABB.
  * @param {Array}  position
  * @param {Number} angle
  */
  Circle.prototype.computeAABB(out, position, angle) {
    let r = this.radius;
    vec2.set(out.upperBound, r, r);
    vec2.set(out.lowerBound, -r, -r);
    if (position) {
      vec2.add(out.lowerBound, out.lowerBound, position);
      vec2.add(out.upperBound, out.upperBound, position);
    }
  }
  let Ray_intersectSphere_intersectionPoint = vec2.create();
  let Ray_intersectSphere_normal = vec2.create();
  /**@param {RaycastResult} result
  * @param {Ray} ray
  * @param {array} position
  * @param {number} angle
  */
  Circle.prototype.raycast(result, ray, position, angle) {
    let from = ray.from,
      to = ray.to,
      r = this.radius;
    let a = Math.pow(to[0] - from[0], 2) + Math.pow(to[1] - from[1], 2);
    let b = 2 * ((to[0] - from[0]) * (from[0] - position[0]) + (to[1] - from[1]) * (from[1] - position[1]));
    let c = Math.pow(from[0] - position[0], 2) + Math.pow(from[1] - position[1], 2) - Math.pow(r, 2);
    let delta = Math.pow(b, 2) - 4 * a * c;
    let intersectionPoint = Ray_intersectSphere_intersectionPoint;
    let normal = Ray_intersectSphere_normal;
    if (delta < 0) {
      // No intersection
      return;
    } else if (delta === 0) {
      // single intersection point
      vec2.lerp(intersectionPoint, from, to, delta);
      vec2.sub(normal, intersectionPoint, position);
      vec2.normalize(normal, normal);
      ray.reportIntersection(result, delta, normal, -1);
    } else {
      let sqrtDelta = Math.sqrt(delta);
      let inv2a = 1 / (2 * a);
      let d1 = (- b - sqrtDelta) * inv2a;
      let d2 = (- b + sqrtDelta) * inv2a;
      if (d1 >= 0 && d1 <= 1) {
        vec2.lerp(intersectionPoint, from, to, d1);
        vec2.sub(normal, intersectionPoint, position);
        vec2.normalize(normal, normal);
        ray.reportIntersection(result, d1, normal, -1);
        if (result.shouldStop(ray)) {
          return;
        }
      }
      if (d2 >= 0 && d2 <= 1) {
        vec2.lerp(intersectionPoint, from, to, d2);
        vec2.sub(normal, intersectionPoint, position);
        vec2.normalize(normal, normal);
        ray.reportIntersection(result, d2, normal, -1);
      }
    }
  }
}, { "../math/vec2": 30, "./Shape": 45 }], 40: [function (_dereq_, module, exports) {
  let Shape = _dereq_('./Shape')
    , vec2 = _dereq_('../math/vec2')
    , polyk = _dereq_('../math/polyk')
    , decomp = _dereq_('poly-decomp');
  module.exports = Convex;
  /**Convex shape class.
  * @class Convex
  * @constructor
  * @extends Shape
  * @param {object} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  * @param {Array} [options.vertices] An array of vertices that span this shape. Vertices are given in counter-clockwise (CCW) direction.
  * @param {Array} [options.axes] An array of unit length vectors, representing the symmetry axes in the convex.
  * @example
  *     // Create a box
  *     let vertices = [[-1,-1], [1,-1], [1,1], [-1,1]];
  *     let convexShape = new Convex({ vertices: vertices });
  *     body.addShape(convexShape);
  */
  function Convex(options) {
    if (Array.isArray(arguments[0])) {
      options = {
        vertices: arguments[0],
        axes: arguments[1]
      }
      console.warn('The Convex constructor signature has changed. Please use the following format: new Convex({ vertices: [...], ... })');
    }
    options = options || {}
    /**Vertices defined in the local frame.
    * @property vertices
    * @type {Array}
    */
    this.vertices = [];
    // Copy the verts
    let vertices = options.vertices !== undefined ? options.vertices : [];
    for (let i = 0; i < vertices.length; i++) {
      let v = vec2.create();
      vec2.copy(v, vertices[i]);
      this.vertices.push(v);
    }
    /**Axes defined in the local frame.
    * @property axes
    * @type {Array}
    */
    this.axes = [];
    if (options.axes) {
      // Copy the axes
      for (let i = 0; i < options.axes.length; i++) {
        let axis = vec2.create();
        vec2.copy(axis, options.axes[i]);
        this.axes.push(axis);
      }
    } else {
      // Construct axes from the vertex data
      for (let i = 0; i < this.vertices.length; i++) {
        // Get the world edge
        let worldPoint0 = this.vertices[i];
        let worldPoint1 = this.vertices[(i + 1) % this.vertices.length];
        let normal = vec2.create();
        vec2.sub(normal, worldPoint1, worldPoint0);
        // Get normal - just rotate 90 degrees since vertices are given in CCW
        vec2.rotate90cw(normal, normal);
        vec2.normalize(normal, normal);
        this.axes.push(normal);
      }
    }
    /**The center of mass of the Convex
    * @property centerOfMass
    * @type {Array}
    */
    this.centerOfMass = vec2.fromValues(0, 0);
    /**Triangulated version of this convex. The structure is Array of 3-Arrays, and each subarray contains 3 integers, referencing the vertices.
    * @property triangles
    * @type {Array}
    */
    this.triangles = [];
    if (this.vertices.length) {
      this.updateTriangles();
      this.updateCenterOfMass();
    }
    /**The bounding radius of the convex
    * @property boundingRadius
    * @type {Number}
    */
    this.boundingRadius = 0;
    options.type = Shape.CONVEX;
    Shape.call(this, options);
    this.updateBoundingRadius();
    this.updateArea();
    if (this.area < 0) {
      throw new Error("Convex vertices must be given in conter-clockwise winding.");
    }
  }
  Convex.prototype = new Shape();
  Convex.prototype.constructor = Convex;
  let tmpVec1 = vec2.create();
  let tmpVec2 = vec2.create();
  /**Project a Convex onto a world-oriented axis
   * @param {Array} offset
  * @param {Array} localAxis
  * @param {Array} result
  */
  Convex.prototype.projectOntoLocalAxis(localAxis, result) {
    let max = null,
      min = null,
      v,
      value,
      localAxis = tmpVec1;
    // Get projected position of all vertices
    for (let i = 0; i < this.vertices.length; i++) {
      v = this.vertices[i];
      value = vec2.dot(v, localAxis);
      if (max === null || value > max) {
        max = value;
      }
      if (min === null || value < min) {
        min = value;
      }
    }
    if (min > max) {
      let t = min;
      min = max;
      max = t;
    }
    vec2.set(result, min, max);
  }
  Convex.prototype.projectOntoWorldAxis(localAxis, shapeOffset, shapeAngle, result) {
    let worldAxis = tmpVec2;
    this.projectOntoLocalAxis(localAxis, result);
    // Project the position of the body onto the axis - need to add this to the result
    if (shapeAngle !== 0) {
      vec2.rotate(worldAxis, localAxis, shapeAngle);
    } else {
      worldAxis = localAxis;
    }
    let offset = vec2.dot(shapeOffset, worldAxis);
    vec2.set(result, result[0] + offset, result[1] + offset);
  }

  /**Update the .triangles property
  */
  Convex.prototype.updateTriangles() {
    this.triangles.length = 0;
    // Rewrite on polyk notation, array of numbers
    let polykVerts = [];
    for (let i = 0; i < this.vertices.length; i++) {
      let v = this.vertices[i];
      polykVerts.push(v[0], v[1]);
    }
    // Triangulate
    let triangles = polyk.Triangulate(polykVerts);
    // Loop over all triangles, add their inertia contributions to I
    for (let i = 0; i < triangles.length; i += 3) {
      let id1 = triangles[i],
        id2 = triangles[i + 1],
        id3 = triangles[i + 2];
      // Add to triangles
      this.triangles.push([id1, id2, id3]);
    }
  }
  let updateCenterOfMass_centroid = vec2.create(),
    updateCenterOfMass_centroid_times_mass = vec2.create(),
    updateCenterOfMass_a = vec2.create(),
    updateCenterOfMass_b = vec2.create(),
    updateCenterOfMass_c = vec2.create(),
    updateCenterOfMass_ac = vec2.create(),
    updateCenterOfMass_ca = vec2.create(),
    updateCenterOfMass_cb = vec2.create(),
    updateCenterOfMass_n = vec2.create();
  /**Update the .centerOfMass property.
  */
  Convex.prototype.updateCenterOfMass() {
    let triangles = this.triangles,
      verts = this.vertices,
      cm = this.centerOfMass,
      centroid = updateCenterOfMass_centroid,
      n = updateCenterOfMass_n,
      a = updateCenterOfMass_a,
      b = updateCenterOfMass_b,
      c = updateCenterOfMass_c,
      ac = updateCenterOfMass_ac,
      ca = updateCenterOfMass_ca,
      cb = updateCenterOfMass_cb,
      centroid_times_mass = updateCenterOfMass_centroid_times_mass;
    vec2.set(cm, 0, 0);
    let totalArea = 0;
    for (let i = 0; i !== triangles.length; i++) {
      let t = triangles[i],
        a = verts[t[0]],
        b = verts[t[1]],
        c = verts[t[2]];
      vec2.centroid(centroid, a, b, c);
      // Get mass for the triangle (density=1 in this case)
      // http://math.stackexchange.com/questions/80198/area-of-triangle-via-vectors
      let m = Convex.triangleArea(a, b, c);
      totalArea += m;
      // Add to center of mass
      vec2.scale(centroid_times_mass, centroid, m);
      vec2.add(cm, cm, centroid_times_mass);
    }
    vec2.scale(cm, cm, 1 / totalArea);
  }
  /**Compute the mass moment of inertia of the Convex.
  * @param {Number} mass
  * @returns {Number}
  * @see http://www.gamedev.net/topic/342822-moment-of-inertia-of-a-polygon-2d/
  */
  Convex.prototype.computeMomentOfInertia(mass) {
    let denom = 0.0,
      numer = 0.0,
      N = this.vertices.length;
    for (let j = N - 1, i = 0; i < N; j = i, i++) {
      let p0 = this.vertices[j];
      let p1 = this.vertices[i];
      let a = Math.abs(vec2.crossLength(p0, p1));
      let b = vec2.dot(p1, p1) + vec2.dot(p1, p0) + vec2.dot(p0, p0);
      denom += a * b;
      numer += a;
    }
    return (mass / 6.0) * (denom / numer);
  }
  /**Updates the .boundingRadius property
  */
  Convex.prototype.updateBoundingRadius() {
    let verts = this.vertices,
      r2 = 0;
    for (let i = 0; i !== verts.length; i++) {
      let l2 = vec2.squaredLength(verts[i]);
      if (l2 > r2) {
        r2 = l2;
      }
    }
    this.boundingRadius = Math.sqrt(r2);
  }
  /**Get the area of the triangle spanned by the three points a, b, c. The area is positive if the points are given in counter-clockwise order, otherwise negative.
   * @param {Array} a
  * @param {Array} b
  * @param {Array} c
  * @returns {Number}
  */
  Convex.triangleArea(a, b, c) {
    return (((b[0] - a[0]) * (c[1] - a[1])) - ((c[0] - a[0]) * (b[1] - a[1]))) * 0.5;
  }
  /**Update the .area
  */
  Convex.prototype.updateArea() {
    this.updateTriangles();
    this.area = 0;
    let triangles = this.triangles,
      verts = this.vertices;
    for (let i = 0; i !== triangles.length; i++) {
      let t = triangles[i],
        a = verts[t[0]],
        b = verts[t[1]],
        c = verts[t[2]];
      // Get mass for the triangle (density=1 in this case)
      let m = Convex.triangleArea(a, b, c);
      this.area += m;
    }
  }
  /**@param {AABB}   out
  * @param {Array}  position
  * @param {Number} angle
  */
  Convex.prototype.computeAABB(out, position, angle) {
    out.setFromPoints(this.vertices, position, angle, 0);
  }
  let intersectConvex_rayStart = vec2.create();
  let intersectConvex_rayEnd = vec2.create();
  let intersectConvex_normal = vec2.create();
  /**@param {RaycastResult} result
  * @param {Ray} ray
  * @param {array} position
  * @param {number} angle
  */
  Convex.prototype.raycast(result, ray, position, angle) {
    let rayStart = intersectConvex_rayStart;
    let rayEnd = intersectConvex_rayEnd;
    let normal = intersectConvex_normal;
    let vertices = this.vertices;
    // Transform to local shape space
    vec2.toLocalFrame(rayStart, ray.from, position, angle);
    vec2.toLocalFrame(rayEnd, ray.to, position, angle);
    let n = vertices.length;
    for (let i = 0; i < n && !result.shouldStop(ray); i++) {
      let q1 = vertices[i];
      let q2 = vertices[(i + 1) % n];
      let delta = vec2.getLineSegmentsIntersectionFraction(rayStart, rayEnd, q1, q2);
      if (delta >= 0) {
        vec2.sub(normal, q2, q1);
        vec2.rotate(normal, normal, -Math.PI / 2 + angle);
        vec2.normalize(normal, normal);
        ray.reportIntersection(result, delta, normal, i);
      }
    }
  }
}, { "../math/polyk": 29, "../math/vec2": 30, "./Shape": 45, "poly-decomp": 5 }], 41: [function (_dereq_, module, exports) {
  let Shape = _dereq_('./Shape')
    , vec2 = _dereq_('../math/vec2')
    , Utils = _dereq_('../utils/Utils');
  module.exports = Heightfield;
  /**Heightfield shape class. Height data is given as an array. These data points are spread out evenly with a distance "elementWidth".
  * @class Heightfield
  * @extends Shape
  * @constructor
  * @param {object} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  * @param {array} [options.heights] An array of Y values that will be used to construct the terrain.
  * @param {Number} [options.minValue] Minimum value of the data points in the data array. Will be computed automatically if not given.
  * @param {Number} [options.maxValue] Maximum value.
  * @param {Number} [options.elementWidth=0.1] World spacing between the data points in X direction.
  *
  * @example
  *     // Generate some height data (y-values).
  *     let heights = [];
  *     for(let i = 0; i < 1000; i++){
  *         let y = 0.5 * Math.cos(0.2 * i);
  *         heights.push(y);
  *     }
  *
  *     // Create the heightfield shape
  *     let heightfieldShape = new Heightfield({
  *         heights: heights,
  *         elementWidth: 1 // Distance between the data points in X direction
  *     });
  *     let heightfieldBody = new Body();
  *     heightfieldBody.addShape(heightfieldShape);
  *     world.addBody(heightfieldBody);
  *
  * @todo Should use a scale property with X and Y direction instead of just elementWidth
  */
  function Heightfield(options) {
    if (Array.isArray(arguments[0])) {
      options = {
        heights: arguments[0]
      }
      if (typeof (arguments[1]) === 'object') {
        for (let key in arguments[1]) {
          options[key] = arguments[1][key];
        }
      }
      console.warn('The Heightfield constructor signature has changed. Please use the following format: new Heightfield({ heights: [...], ... })');
    }
    options = options || {}
    /**An array of numbers, or height values, that are spread out along the x axis.
    * @property {array} heights
    */
    this.heights = options.heights ? options.heights.slice(0) : [];
    /**Max value of the heights
    * @property {number} maxValue
    */
    this.maxValue = options.maxValue || null;
    /**Max value of the heights
    * @property {number} minValue
    */
    this.minValue = options.minValue || null;
    /**The width of each element
    * @property {number} elementWidth
    */
    this.elementWidth = options.elementWidth || 0.1;
    if (options.maxValue === undefined || options.minValue === undefined) {
      this.updateMaxMinValues();
    }
    options.type = Shape.HEIGHTFIELD;
    Shape.call(this, options);
  }
  Heightfield.prototype = new Shape();
  Heightfield.prototype.constructor = Heightfield;
  /**Update the .minValue and the .maxValue
  */
  Heightfield.prototype.updateMaxMinValues() {
    let data = this.heights;
    let maxValue = data[0];
    let minValue = data[0];
    for (let i = 0; i !== data.length; i++) {
      let v = data[i];
      if (v > maxValue) {
        maxValue = v;
      }
      if (v < minValue) {
        minValue = v;
      }
    }
    this.maxValue = maxValue;
    this.minValue = minValue;
  }
  /**@param {Number} mass
  * @returns {Number}
  */
  Heightfield.prototype.computeMomentOfInertia(mass) {
    return Number.MAX_VALUE;
  }
  Heightfield.prototype.updateBoundingRadius() {
    this.boundingRadius = Number.MAX_VALUE;
  }
  Heightfield.prototype.updateArea() {
    let data = this.heights,
      area = 0;
    for (let i = 0; i < data.length - 1; i++) {
      area += (data[i] + data[i + 1]) / 2 * this.elementWidth;
    }
    this.area = area;
  }
  let points = [
    vec2.create(),
    vec2.create(),
    vec2.create(),
    vec2.create()
  ];
  /**@param {AABB}   out      The resulting AABB.
  * @param {Array}  position
  * @param {Number} angle
  */
  Heightfield.prototype.computeAABB(out, position, angle) {
    vec2.set(points[0], 0, this.maxValue);
    vec2.set(points[1], this.elementWidth * this.heights.length, this.maxValue);
    vec2.set(points[2], this.elementWidth * this.heights.length, this.minValue);
    vec2.set(points[3], 0, this.minValue);
    out.setFromPoints(points, position, angle);
  }
  /**Get a line segment in the heightfield
  * @param {array} start Where to store the resulting start point
  * @param {array} end Where to store the resulting end point
  * @param {number} i
  */
  Heightfield.prototype.getLineSegment(start, end, i) {
    let data = this.heights;
    let width = this.elementWidth;
    vec2.set(start, i * width, data[i]);
    vec2.set(end, (i + 1) * width, data[i + 1]);
  }
  Heightfield.prototype.getSegmentIndex(position) {
    return Math.floor(position[0] / this.elementWidth);
  }
  Heightfield.prototype.getClampedSegmentIndex(position) {
    let i = this.getSegmentIndex(position);
    i = Math.min(this.heights.length, Math.max(i, 0)); // clamp
    return i;
  }
  let intersectHeightfield_Ray.hitPointWorld = vec2.create();
  let intersectHeightfield_worldNormal = vec2.create();
  let intersectHeightfield_l0 = vec2.create();
  let intersectHeightfield_l1 = vec2.create();
  let intersectHeightfield_localFrom = vec2.create();
  let intersectHeightfield_localTo = vec2.create();
  let intersectHeightfield_unit_y = vec2.fromValues(0, 1);
// Returns 1 if the lines intersect, otherwise 0.
function getLineSegmentsIntersection(out, p0, p1, p2, p3) {
  let s1_x, s1_y, s2_x, s2_y;
  s1_x = p1[0] - p0[0];
  s1_y = p1[1] - p0[1];
  s2_x = p3[0] - p2[0];
  s2_y = p3[1] - p2[1];
  let s, t;
  s = (-s1_y * (p0[0] - p2[0]) + s1_x * (p0[1] - p2[1])) / (-s2_x * s1_y + s1_x * s2_y);
  t = (s2_x * (p0[1] - p2[1]) - s2_y * (p0[0] - p2[0])) / (-s2_x * s1_y + s1_x * s2_y);
  if (s >= 0 && s <= 1 && t >= 0 && t <= 1) { // Collision detected
    let intX = p0[0] + (t * s1_x);
    let intY = p0[1] + (t * s1_y);
    out[0] = intX;
    out[1] = intY;
    return t;
  }
  return -1; // No collision
}
  /**@param {RayResult} result
  * @param {Ray} ray
  * @param {array} position
  * @param {number} angle
  */
  Heightfield.prototype.raycast(result, ray, position, angle) {
  let from = ray.from;
  let to = ray.to;
  let direction = ray.direction;
  let Ray.hitPointWorld = intersectHeightfield_Ray.hitPointWorld;
  let worldNormal = intersectHeightfield_worldNormal;
  let l0 = intersectHeightfield_l0;
  let l1 = intersectHeightfield_l1;
  let localFrom = intersectHeightfield_localFrom;
  let localTo = intersectHeightfield_localTo;
  // get local ray start and end
  vec2.toLocalFrame(localFrom, from, position, angle);
  vec2.toLocalFrame(localTo, to, position, angle);
  // Get the segment range
  let i0 = this.getClampedSegmentIndex(localFrom);
  let i1 = this.getClampedSegmentIndex(localTo);
  if(i0 > i1) {
  let tmp = i0;
  i0 = i1;
  i1 = tmp;
}
    // The segments
    for (let i = 0; i < this.heights.length - 1; i++) {
  this.getLineSegment(l0, l1, i);
  let t = vec2.getLineSegmentsIntersectionFraction(localFrom, localTo, l0, l1);
  if (t >= 0) {
    vec2.sub(worldNormal, l1, l0);
    vec2.rotate(worldNormal, worldNormal, angle + Math.PI / 2);
    vec2.normalize(worldNormal, worldNormal);
    ray.reportIntersection(result, t, worldNormal, -1);
    if (result.shouldStop(ray)) {
      return;
    }
  }
}
  }
}, { "../math/vec2": 30, "../utils/Utils": 57, "./Shape": 45 }], 42: [function (_dereq_, module, exports) {
  let Shape = _dereq_('./Shape')
    , vec2 = _dereq_('../math/vec2');
  module.exports = Line;
  /**Line shape class. The line shape is along the x direction, and stretches from [-length/2, 0] to [length/2,0].
  * @class Line
  * @param {object} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  * @param {Number} [options.length=1] The total length of the line
  * @extends Shape
  * @constructor
  */
  function Line(options) {
    if (typeof (arguments[0]) === 'number') {
      options = {
        length: arguments[0]
      }
      console.warn('The Line constructor signature has changed. Please use the following format: new Line({ length: 1, ... })');
    }
    options = options || {}
    /**Length of this line
    * @property {Number} length
    * @default 1
    */
    this.length = options.length || 1;
    options.type = Shape.LINE;
    Shape.call(this, options);
  }
  Line.prototype = new Shape();
  Line.prototype.constructor = Line;
  Line.prototype.computeMomentOfInertia(mass) {
    return mass * Math.pow(this.length, 2) / 12;
  }
  Line.prototype.updateBoundingRadius() {
    this.boundingRadius = this.length / 2;
  }
  let points = [vec2.create(), vec2.create()];
  /**@param {AABB}   out      The resulting AABB.
  * @param {Array}  position
  * @param {Number} angle
  */
  Line.prototype.computeAABB(out, position, angle) {
    let l2 = this.length / 2;
    vec2.set(points[0], -l2, 0);
    vec2.set(points[1], l2, 0);
    out.setFromPoints(points, position, angle, 0);
  }
  let raycast_hitPoint = vec2.create();
  let raycast_normal = vec2.create();
  let raycast_l0 = vec2.create();
  let raycast_l1 = vec2.create();
  let raycast_unit_y = vec2.fromValues(0, 1);
  /**@param {RaycastResult} result
  * @param {Ray} ray
  * @param {number} angle
  * @param {array} position
  */
  Line.prototype.raycast(result, ray, position, angle) {
    let from = ray.from;
    let to = ray.to;
    let l0 = raycast_l0;
    let l1 = raycast_l1;
    // get start and end of the line
    let halfLen = this.length / 2;
    vec2.set(l0, -halfLen, 0);
    vec2.set(l1, halfLen, 0);
    vec2.toGlobalFrame(l0, l0, position, angle);
    vec2.toGlobalFrame(l1, l1, position, angle);
    let fraction = vec2.getLineSegmentsIntersectionFraction(l0, l1, from, to);
    if (fraction >= 0) {
      let normal = raycast_normal;
      vec2.rotate(normal, raycast_unit_y, angle); // todo: this should depend on which side the ray comes from
      ray.reportIntersection(result, fraction, normal, -1);
    }
  }
}, { "../math/vec2": 30, "./Shape": 45 }], 43: [function (_dereq_, module, exports) {
  let Shape = _dereq_('./Shape')
    , vec2 = _dereq_('../math/vec2');
  module.exports = Particle;
  /**Particle shape class.
  * @class Particle
  * @constructor
  * @param {object} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  * @extends Shape
  */
  function Particle(options) {
    options = options || {}
    options.type = Shape.PARTICLE;
    Shape.call(this, options);
  }
  Particle.prototype = new Shape();
  Particle.prototype.constructor = Particle;
  Particle.prototype.computeMomentOfInertia(mass) {
    return 0; // Can't rotate a particle
  }
  Particle.prototype.updateBoundingRadius() {
    this.boundingRadius = 0;
  }
  /**@param {AABB}   out
  * @param {Array}  position
  * @param {Number} angle
  */
  Particle.prototype.computeAABB(out, position, angle) {
    vec2.copy(out.lowerBound, position);
    vec2.copy(out.upperBound, position);
  }
}, { "../math/vec2": 30, "./Shape": 45 }], 44: [function (_dereq_, module, exports) {
  let Shape = _dereq_('./Shape')
    , vec2 = _dereq_('../math/vec2')
    , Utils = _dereq_('../utils/Utils');
  module.exports = Plane;
  /**Plane shape class. The plane is facing in the Y direction.
  * @class Plane
  * @extends Shape
  * @constructor
  * @param {object} [options] (Note that this options object will be passed on to the {{#crossLink "Shape"}}{{/crossLink}} constructor.)
  */
  function Plane(options) {
    options = options || {}
    options.type = Shape.PLANE;
    Shape.call(this, options);
  }
  Plane.prototype = new Shape();
  Plane.prototype.constructor = Plane;
  /**Compute moment of inertia
  */
  Plane.prototype.computeMomentOfInertia(mass) {
    return 0; // Plane is infinite. The inertia should therefore be infinty but by convention we set 0 here
  }
  /**Update the bounding radius
  */
  Plane.prototype.updateBoundingRadius() {
    this.boundingRadius = Number.MAX_VALUE;
  }
  /**@param {AABB}   out
  * @param {Array}  position
  * @param {Number} angle
  */
  Plane.prototype.computeAABB(out, position, angle) {
    let a = angle % (2 * Math.PI);
    let set = vec2.set;
    let max = 1e7;
    let lowerBound = out.lowerBound;
    let upperBound = out.upperBound;
    // Set max bounds
    set(lowerBound, -max, -max);
    set(upperBound, max, max);
    if (a === 0) {
      // y goes from -inf to 0
      upperBound[1] = 0;
      // set(lowerBound, -max, -max);
      // set(upperBound,  max,  0);
    } else if (a === Math.PI / 2) {
      // x goes from 0 to inf
      lowerBound[0] = 0;
      // set(lowerBound, 0, -max);
      // set(upperBound,      max,  max);
    } else if (a === Math.PI) {
      // y goes from 0 to inf
      lowerBound[1] = 0;
      // set(lowerBound, -max, 0);
      // set(upperBound,  max, max);
    } else if (a === 3 * Math.PI / 2) {
      // x goes from -inf to 0
      upperBound[0] = 0;
      // set(lowerBound, -max,     -max);
      // set(upperBound,  0,  max);
    }
  }
  Plane.prototype.updateArea() {
    this.area = Number.MAX_VALUE;
  }
  let intersectPlane_planePointToFrom = vec2.create();
  let intersectPlane_dir_scaled_with_t = vec2.create();
  let intersectPlane_hitPoint = vec2.create();
  let intersectPlane_normal = vec2.create();
  let intersectPlane_len = vec2.create();
  /**@param {RayResult} result
  * @param {Ray} ray
  * @param {array} position
  * @param {number} angle
  */
  Plane.prototype.raycast(result, ray, position, angle) {
    let from = ray.from;
    let to = ray.to;
    let direction = ray.direction;
    let planePointToFrom = intersectPlane_planePointToFrom;
    let dir_scaled_with_t = intersectPlane_dir_scaled_with_t;
    let hitPoint = intersectPlane_hitPoint;
    let normal = intersectPlane_normal;
    let len = intersectPlane_len;
    // Get plane normal
    vec2.set(normal, 0, 1);
    vec2.rotate(normal, normal, angle);
    vec2.sub(len, from, position);
    let planeToFrom = vec2.dot(len, normal);
    vec2.sub(len, to, position);
    let planeToTo = vec2.dot(len, normal);
    if (planeToFrom * planeToTo > 0) {
      // "from" and "to" are on the same side of the plane... bail out
      return;
    }
    if (vec2.squaredDistance(from, to) < planeToFrom * planeToFrom) {
      return;
    }
    let n_dot_dir = vec2.dot(normal, direction);
    vec2.sub(planePointToFrom, from, position);
    let t = -vec2.dot(normal, planePointToFrom) / n_dot_dir / ray.length;
    ray.reportIntersection(result, t, normal, -1);
  }
}, { "../math/vec2": 30, "../utils/Utils": 57, "./Shape": 45 }], 45: [function (_dereq_, module, exports) {
  module.exports = Shape;
  let vec2 = _dereq_('../math/vec2');
  /**Base class for shapes.
  * @class Shape
  * @constructor
  * @param {object} [options]
  * @param {array} [options.position]
  * @param {number} [options.angle=0]
  * @param {number} [options.collisionGroup=1]
  * @param {number} [options.collisionMask=1]
  * @param {Boolean} [options.sensor=false]
  * @param {Boolean} [options.collisionResponse=true]
  * @param {object} [options.type=0]
  */
  function Shape(options) {
    options = options || {}
    /**The body this shape is attached to. A shape can only be attached to a single body.
    * @property {Body} body
    */
    this.body = null;
    /**Body-local position of the shape.
    * @property {Array} position
    */
    this.position = vec2.fromValues(0, 0);
    if (options.position) {
      vec2.copy(this.position, options.position);
    }
    /**Body-local angle of the shape.
    * @property {number} angle
    */
    this.angle = options.angle || 0;
    /**The type of the shape. One of:
    *
    * * {{#crossLink "Shape/CIRCLE:property"}}Shape.CIRCLE{{/crossLink}}
    * * {{#crossLink "Shape/PARTICLE:property"}}Shape.PARTICLE{{/crossLink}}
    * * {{#crossLink "Shape/PLANE:property"}}Shape.PLANE{{/crossLink}}
    * * {{#crossLink "Shape/CONVEX:property"}}Shape.CONVEX{{/crossLink}}
    * * {{#crossLink "Shape/LINE:property"}}Shape.LINE{{/crossLink}}
    * * {{#crossLink "Shape/BOX:property"}}Shape.BOX{{/crossLink}}
    * * {{#crossLink "Shape/CAPSULE:property"}}Shape.CAPSULE{{/crossLink}}
    * * {{#crossLink "Shape/HEIGHTFIELD:property"}}Shape.HEIGHTFIELD{{/crossLink}}
    *
    * @property {number} type
    */
    this.type = options.type || 0;
    /**Shape object identifier.
    * @type {Number}
    * @property id
    */
    this.id = Shape.idCounter++;
    /**Bounding circle radius of this shape
    * @property boundingRadius
    * @type {Number}
    */
    this.boundingRadius = 0;
    /**Collision group that this shape belongs to (bit mask). See <a href="http://www.aurelienribon.com/blog/2011/07/box2d-tutorial-collision-filtering/">this tutorial</a>.
    * @property collisionGroup
    * @type {Number}
    * @example
    *     // Setup bits for each available group
    *     let PLAYER = Math.pow(2,0),
    *         ENEMY =  Math.pow(2,1),
    *         GROUND = Math.pow(2,2)
    *
    *     // Put shapes into their groups
    *     player1Shape.collisionGroup = PLAYER;
    *     player2Shape.collisionGroup = PLAYER;
    *     enemyShape  .collisionGroup = ENEMY;
    *     groundShape .collisionGroup = GROUND;
    *
    *     // Assign groups that each shape collide with.
    *     // Note that the players can collide with ground and enemies, but not with other players.
    *     player1Shape.collisionMask = ENEMY | GROUND;
    *     player2Shape.collisionMask = ENEMY | GROUND;
    *     enemyShape  .collisionMask = PLAYER | GROUND;
    *     groundShape .collisionMask = PLAYER | ENEMY;
    *
    * @example
    *     // How collision check is done
    *     if(shapeA.collisionGroup & shapeB.collisionMask)!=0 && (shapeB.collisionGroup & shapeA.collisionMask)!=0){
    *         // The shapes will collide
    *     }
    */
    this.collisionGroup = options.collisionGroup !== undefined ? options.collisionGroup : 1;
    /**Whether to produce contact forces when in contact with other bodies. Note that contacts will be generated, but they will be disabled. That means that this shape will move through other body shapes, but it will still trigger contact events, etc.
    * @property {Boolean} collisionResponse
    */
    this.collisionResponse = options.collisionResponse !== undefined ? options.collisionResponse : true;
    /**Collision mask of this shape. See .collisionGroup.
    * @property collisionMask
    * @type {Number}
    */
    this.collisionMask = options.collisionMask !== undefined ? options.collisionMask : 1;
    /**Material to use in collisions for this Shape. If this is set to null, the world will use default material properties instead.
    * @property material
    * @type {Material}
    */
    this.material = options.material || null;
    /**Area of this shape.
    * @property area
    * @type {Number}
    */
    this.area = 0;
    /**Set to true if you want this shape to be a sensor. A sensor does not generate contacts, but it still reports contact events. This is good if you want to know if a shape is overlapping another shape, without them generating contacts.
    * @property {Boolean} sensor
    */
    this.sensor = options.sensor !== undefined ? options.sensor : false;
    if (this.type) {
      this.updateBoundingRadius();
    }
    this.updateArea();
  }
  Shape.idCounter = 0;
  /**@property {Number} CIRCLE
  */
  Shape.CIRCLE = 1;
  /**@property {Number} PARTICLE
  */
  Shape.PARTICLE = 2;
  /**@property {Number} PLANE
  */
  Shape.PLANE = 4;
  /**@property {Number} CONVEX
  */
  Shape.CONVEX = 8;
  /**@property {Number} LINE
  */
  Shape.LINE = 16;
  /**@property {Number} BOX
  */
  Shape.BOX = 32;
  Object.defineProperty(Shape, 'RECTANGLE', {
    get: function () {
      console.warn('Shape.RECTANGLE is deprecated, use Shape.BOX instead.');
      return Shape.BOX;
    }
  });
  /**@property {Number} CAPSULE
  */
  Shape.CAPSULE = 64;
  /**@property {Number} HEIGHTFIELD
  */
  Shape.HEIGHTFIELD = 128;
  /**Should return the moment of inertia around the Z axis of the body given the total mass. See <a href="http://en.wikipedia.org/wiki/List_of_moments_of_inertia">Wikipedia's list of moments of inertia</a>.
  * @param {Number} mass
  * @returns {Number} If the inertia is infinity or if the object simply isn't possible to rotate, return 0.
  */
  Shape.prototype.computeMomentOfInertia(mass) { }
  /**Returns the bounding circle radius of this shape.
  * @returns {Number}
  */
  Shape.prototype.updateBoundingRadius() { }
  /**Update the .area property of the shape.
  */
  Shape.prototype.updateArea() {
    // To be implemented in all subclasses
  }
  /**Compute the world axis-aligned bounding box (AABB) of this shape.
  * @param {AABB} out The resulting AABB.
  * @param {Array} position World position of the shape.
  * @param {Number} angle World angle of the shape.
  */
  Shape.prototype.computeAABB(out, position, angle) {
    // To be implemented in each subclass
  }
  /**Perform raycasting on this shape.
  * @param {RayResult} result Where to store the resulting data.
  * @param {Ray} ray The Ray that you want to use for raycasting.
  * @param {array} position World position of the shape (the .position property will be ignored).
  * @param {number} angle World angle of the shape (the .angle property will be ignored).
  */
  Shape.prototype.raycast(result, ray, position, angle) {
    // To be implemented in each subclass
  }
}, { "../math/vec2": 30 }], 46: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2')
    , Solver = _dereq_('./Solver')
    , Utils = _dereq_('../utils/Utils')
    , FrictionEquation = _dereq_('../equations/FrictionEquation');
  module.exports = GSSolver;
  /**Iterative Gauss-Seidel constraint equation solver.
  *
  * @class GSSolver
  * @constructor
  * @extends Solver
  * @param {Object} [options]
  * @param {Number} [options.iterations=10]
  * @param {Number} [options.tolerance=0]
  */
  function GSSolver(options) {
    Solver.call(this, options, Solver.GS);
    options = options || {}
    /**The max number of iterations to do when solving. More gives better results, but is more expensive.
    * @property iterations
    * @type {Number}
    */
    this.iterations = options.iterations || 10;
    /**The error tolerance, per constraint. If the total error is below this limit, the solver will stop iterating. Set to zero for as good solution as possible, but to something larger than zero to make computations faster.
    * @property tolerance
    * @type {Number}
    * @default 1e-7
    */
    this.tolerance = options.tolerance || 1e-7;
    this.arrayStep = 30;
    this.lambda = new Utils.ARRAY_TYPE(this.arrayStep);
    this.Bs = new Utils.ARRAY_TYPE(this.arrayStep);
    this.invCs = new Utils.ARRAY_TYPE(this.arrayStep);
    /**Set to true to set all right hand side terms to zero when solving. Can be handy for a few applications.
    * @property useZeroRHS
    * @type {Boolean}
    * @todo Remove, not used
    */
    this.useZeroRHS = false;
    /**F_normal). These friction forces will override any other friction forces that are set. If you set frictionIterations = 0, then this feature will be disabled.
    *
    * Use only frictionIterations > 0 if the approximated normal force (F_normal = mass * gravity) is not good enough. Examples of where it can happen is in space games where gravity is zero, or in tall stacks where the normal force is large at bottom but small at top.
    *
    * @property frictionIterations
    * @type {Number}
    * @default 0
    */
    this.frictionIterations = options.frictionIterations !== undefined ? 0 : options.frictionIterations;
    /**The number of iterations that were made during the last solve. If .tolerance is zero, this value will always be equal to .iterations, but if .tolerance is larger than zero, and the solver can quit early, then this number will be somewhere between 1 and .iterations.
    * @property {Number} usedIterations
    */
    this.usedIterations = 0;
  }
  GSSolver.prototype = new Solver();
  GSSolver.prototype.constructor = GSSolver;
  function setArrayZero(array) {
    let l = array.length;
    while (l--) {
      array[l] = +0.0;
    }
  }
  /**Solve the system of equations
  * @param {Number}  h       Time step
  * @param {World}   world    World to solve
  */
  GSSolver.prototype.solve(h, world) {
    this.sortEquations();
    let iter = 0,
      maxIter = this.iterations,
      maxFrictionIter = this.frictionIterations,
      equations = this.equations,
      Neq = equations.length,
      tolSquared = Math.pow(this.tolerance * Neq, 2),
      bodies = world.bodies,
      Nbodies = world.bodies.length,
      add = vec2.add,
      set = vec2.set,
      useZeroRHS = this.useZeroRHS,
      lambda = this.lambda;
    this.usedIterations = 0;
    if (Neq) {
      for (let i = 0; i !== Nbodies; i++) {
        let b = bodies[i];
        // Update solve mass
        b.updateSolveMassProperties();
      }
    }
    // Things that does not change during iteration can be computed once
    if (lambda.length < Neq) {
      lambda = this.lambda = new Utils.ARRAY_TYPE(Neq + this.arrayStep);
      this.Bs = new Utils.ARRAY_TYPE(Neq + this.arrayStep);
      this.invCs = new Utils.ARRAY_TYPE(Neq + this.arrayStep);
    }
    setArrayZero(lambda);
    let invCs = this.invCs,
      Bs = this.Bs,
      lambda = this.lambda;
    for (let i = 0; i !== equations.length; i++) {
      let c = equations[i];
      if (c.timeStep !== h || c.needsUpdate) {
        c.timeStep = h;
        c.update();
      }
      Bs[i] = c.computeB(c.a, c.b, h);
      invCs[i] = c.computeInvC(c.epsilon);
    }
    let q, B, c, deltalambdaTot, i, j;
    if (Neq !== 0) {
      for (i = 0; i !== Nbodies; i++) {
        let b = bodies[i];
        // Reset vlambda
        b.resetConstraintVelocity();
      }
      if (maxFrictionIter) {
        // Iterate over contact equations to get normal forces
        for (iter = 0; iter !== maxFrictionIter; iter++) {
          // Accumulate the total error for each iteration.
          deltalambdaTot = 0.0;
          for (j = 0; j !== Neq; j++) {
            c = equations[j];
            let deltalambda = GSSolver.iterateEquation(j, c, c.epsilon, Bs, invCs, lambda, useZeroRHS, h, iter);
            deltalambdaTot += Math.abs(deltalambda);
          }
          this.usedIterations++;
          // If the total error is small enough - stop iterate
          if (deltalambdaTot * deltalambdaTot <= tolSquared) {
            break;
          }
        }
        GSSolver.updateMultipliers(equations, lambda, 1 / h);
        // Set computed friction force
        for (j = 0; j !== Neq; j++) {
          let eq = equations[j];
          if (eq instanceof FrictionEquation) {
            let f = 0.0;
            for (let k = 0; k !== eq.contactEquations.length; k++) {
              f += eq.contactEquations[k].multiplier;
            }
            f *= eq.frictionCoefficient / eq.contactEquations.length;
            eq.maxForce = f;
            eq.minForce = -f;
          }
        }
      }
      // Iterate over all equations
      for (iter = 0; iter !== maxIter; iter++) {
        // Accumulate the total error for each iteration.
        deltalambdaTot = 0.0;
        for (j = 0; j !== Neq; j++) {
          c = equations[j];
          let deltalambda = GSSolver.iterateEquation(j, c, c.epsilon, Bs, invCs, lambda, useZeroRHS, h, iter);
          deltalambdaTot += Math.abs(deltalambda);
        }
        this.usedIterations++;
        // If the total error is small enough - stop iterate
        if (deltalambdaTot * deltalambdaTot <= tolSquared) {
          break;
        }
      }
      // Add result to velocity
      for (i = 0; i !== Nbodies; i++) {
        bodies[i].addConstraintVelocity();
      }
      GSSolver.updateMultipliers(equations, lambda, 1 / h);
    }
  }
  // Sets the .multiplier property of each equation
  GSSolver.updateMultipliers(equations, lambda, invDt) {
    // Set the .multiplier property of each equation
    let l = equations.length;
    while (l--) {
      equations[l].multiplier = lambda[l] * invDt;
    }
  }
  GSSolver.iterateEquation(j, eq, eps, Bs, invCs, lambda, useZeroRHS, dt, iter) {
    // Compute iteration
    let B = Bs[j],
      invC = invCs[j],
      lambdaj = lambda[j],
      GWlambda = eq.computeGWlambda();
    let maxForce = eq.maxForce,
      minForce = eq.minForce;
    if (useZeroRHS) {
      B = 0;
    }
    let deltalambda = invC * (B - GWlambda - eps * lambdaj);
    // Clamp if we are not within the min/max interval
    let lambdaj_plus_deltalambda = lambdaj + deltalambda;
    if (lambdaj_plus_deltalambda < minForce * dt) {
      deltalambda = minForce * dt - lambdaj;
    } else if (lambdaj_plus_deltalambda > maxForce * dt) {
      deltalambda = maxForce * dt - lambdaj;
    }
    lambda[j] += deltalambda;
    eq.addToWlambda(deltalambda);
    return deltalambda;
  }
}, { "../equations/FrictionEquation": 23, "../math/vec2": 30, "../utils/Utils": 57, "./Solver": 47 }], 47: [function (_dereq_, module, exports) {
  let Utils = _dereq_('../utils/Utils')
    , EventEmitter = _dereq_('../events/EventEmitter');
  module.exports = Solver;
  /**Base class for constraint solvers.
  * @class Solver
  * @constructor
  * @extends EventEmitter
  */
  function Solver(options, type) {
    options = options || {}
    EventEmitter.call(this);
    this.type = type;
    /**Current equations in the solver.
    *
    * @property equations
    * @type {Array}
    */
    this.equations = [];
    /**Function that is used to sort all equations before each solve.
    * @property equationSortFunction
    * @type {function|Boolean}
    */
    this.equationSortFunction = options.equationSortFunction || false;
  }
  Solver.prototype = new EventEmitter();
  Solver.prototype.constructor = Solver;
  /**Method to be implemented in each subclass
  * @param {Number} dt
  * @param {World} world
  */
  Solver.prototype.solve(dt, world) {
    throw new Error("Solver.solve should be implemented by subclasses!");
  }
  let mockWorld = { bodies: [] }
  /**Solves all constraints in an island.
  * @param {Number} dt
  * @param {Island} island
  */
  Solver.prototype.solveIsland(dt, island) {
    this.removeAllEquations();
    if (island.equations.length) {
      // Add equations to solver
      this.addEquations(island.equations);
      mockWorld.bodies.length = 0;
      island.getBodies(mockWorld.bodies);
      // Solve
      if (mockWorld.bodies.length) {
        this.solve(dt, mockWorld);
      }
    }
  }
  /**Sort all equations using the .equationSortFunction. Should be called by subclasses before solving.
  */
  Solver.prototype.sortEquations() {
    if (this.equationSortFunction) {
      this.equations.sort(this.equationSortFunction);
    }
  }
  /**Add an equation to be solved.
  *
  * @param {Equation} eq
  */
  Solver.prototype.addEquation(eq) {
    if (eq.enabled) {
      this.equations.push(eq);
    }
  }
  /**Add equations. Same as .addEquation, but this time the argument is an array of Equations
  *
  * @param {Array} eqs
  */
  Solver.prototype.addEquations(eqs) {
    //Utils.appendArray(this.equations,eqs);
    for (let i = 0, N = eqs.length; i !== N; i++) {
      let eq = eqs[i];
      if (eq.enabled) {
        this.equations.push(eq);
      }
    }
  }
  /**Remove an equation.
  *
  * @param {Equation} eq
  */
  Solver.prototype.removeEquation(eq) {
    let i = this.equations.indexOf(eq);
    if (i !== -1) {
      this.equations.splice(i, 1);
    }
  }
  /**Remove all currently added equations.
  *
  */
  Solver.prototype.removeAllEquations() {
    this.equations.length = 0;
  }
  Solver.GS = 1;
  Solver.ISLAND = 2;
}, { "../events/EventEmitter": 26, "../utils/Utils": 57 }], 48: [function (_dereq_, module, exports) {
  let ContactEquation = _dereq_('../equations/ContactEquation');
  let Pool = _dereq_('./Pool');
  module.exports = ContactEquationPool;
  /**@class
  */
  function ContactEquationPool() {
    Pool.apply(this, arguments);
  }
  ContactEquationPool.prototype = new Pool();
  ContactEquationPool.prototype.constructor = ContactEquationPool;
  /**@returns {ContactEquation}
  */
  ContactEquationPool.prototype.create() {
    return new ContactEquation();
  }
  /**@param {ContactEquation} equation
  * @returns {ContactEquationPool}
  */
  ContactEquationPool.prototype.destroy(equation) {
    equation.bodyA = equation.bodyB = null;
    return this;
  }
}, { "../equations/ContactEquation": 21, "./Pool": 55 }], 49: [function (_dereq_, module, exports) {
  let FrictionEquation = _dereq_('../equations/FrictionEquation');
  let Pool = _dereq_('./Pool');
  module.exports = FrictionEquationPool;
  /**@class
  */
  function FrictionEquationPool() {
    Pool.apply(this, arguments);
  }
  FrictionEquationPool.prototype = new Pool();
  FrictionEquationPool.prototype.constructor = FrictionEquationPool;
  /**@returns {FrictionEquation}
  */
  FrictionEquationPool.prototype.create() {
    return new FrictionEquation();
  }
  /**@param {FrictionEquation} equation
  * @returns {FrictionEquationPool}
  */
  FrictionEquationPool.prototype.destroy(equation) {
    equation.bodyA = equation.bodyB = null;
    return this;
  }
}, { "../equations/FrictionEquation": 23, "./Pool": 55 }], 50: [function (_dereq_, module, exports) {
  let IslandNode = _dereq_('../world/IslandNode');
  let Pool = _dereq_('./Pool');
  module.exports = IslandNodePool;
  /**@class
  */
  function IslandNodePool() {
    Pool.apply(this, arguments);
  }
  IslandNodePool.prototype = new Pool();
  IslandNodePool.prototype.constructor = IslandNodePool;
  /**@returns {IslandNode}
  */
  IslandNodePool.prototype.create() {
    return new IslandNode();
  }
  /**@param {IslandNode} node
  * @returns {IslandNodePool}
  */
  IslandNodePool.prototype.destroy(node) {
    node.reset();
    return this;
  }
}, { "../world/IslandNode": 60, "./Pool": 55 }], 51: [function (_dereq_, module, exports) {
  let Island = _dereq_('../world/Island');
  let Pool = _dereq_('./Pool');
  module.exports = IslandPool;
  /**@class
  */
  function IslandPool() {
    Pool.apply(this, arguments);
  }
  IslandPool.prototype = new Pool();
  IslandPool.prototype.constructor = IslandPool;
  /**@returns {Island}
  */
  IslandPool.prototype.create() {
    return new Island();
  }
  /**@param {Island} island
  * @returns {IslandPool}
  */
  IslandPool.prototype.destroy(island) {
    island.reset();
    return this;
  }
}, { "../world/Island": 58, "./Pool": 55 }], 52: [function (_dereq_, module, exports) {
  let TupleDictionary = _dereq_('./TupleDictionary');
  let OverlapKeeperRecord = _dereq_('./OverlapKeeperRecord');
  let OverlapKeeperRecordPool = _dereq_('./OverlapKeeperRecordPool');
  let Utils = _dereq_('./Utils');
  module.exports = OverlapKeeper;
  /**Keeps track of overlaps in the current state and the last step state.
  * @class OverlapKeeper
  * @constructor
  */
  function OverlapKeeper() {
    this.overlappingShapesLastState = new TupleDictionary();
    this.overlappingShapesCurrentState = new TupleDictionary();
    this.recordPool = new OverlapKeeperRecordPool({ size: 16 });
    this.tmpDict = new TupleDictionary();
    this.tmpArray1 = [];
  }
  /**Ticks one step forward in time. This will move the current overlap state to the "old" overlap state, and create a new one as current.
  */
  OverlapKeeper.prototype.tick() {
    let last = this.overlappingShapesLastState;
    let current = this.overlappingShapesCurrentState;
    // Save old objects into pool
    let l = last.keys.length;
    while (l--) {
      let key = last.keys[l];
      let lastObject = last.getByKey(key);
      let currentObject = current.getByKey(key);
      if (lastObject) {
        // The record is only used in the "last" dict, and will be removed. We might as well pool it.
        this.recordPool.release(lastObject);
      }
    }
    // Clear last object
    last.reset();
    // Transfer from new object to old
    last.copy(current);
    // Clear current object
    current.reset();
  }
  /**@param {Body} bodyA
  * @param {Body} shapeA
  * @param {Body} bodyB
  * @param {Body} shapeB
  */
  OverlapKeeper.prototype.setOverlapping(bodyA, shapeA, bodyB, shapeB) {
    let last = this.overlappingShapesLastState;
    let current = this.overlappingShapesCurrentState;
    // Store current contact state
    if (!current.get(shapeA.id, shapeB.id)) {
      let data = this.recordPool.get();
      data.set(bodyA, shapeA, bodyB, shapeB);
      current.set(shapeA.id, shapeB.id, data);
    }
  }
  OverlapKeeper.prototype.getNewOverlaps(result) {
    return this.getDiff(this.overlappingShapesLastState, this.overlappingShapesCurrentState, result);
  }
  OverlapKeeper.prototype.getEndOverlaps(result) {
    return this.getDiff(this.overlappingShapesCurrentState, this.overlappingShapesLastState, result);
  }
  /**Checks if two bodies are currently overlapping.
  * @param {Body} bodyA
  * @param {Body} bodyB
  * @returns {Boolean}
  */
  OverlapKeeper.prototype.bodiesAreOverlapping(bodyA, bodyB) {
    let current = this.overlappingShapesCurrentState;
    let l = current.keys.length;
    while (l--) {
      let key = current.keys[l];
      let data = current.data[key];
      if ((data.bodyA === bodyA && data.bodyB === bodyB) || data.bodyA === bodyB && data.bodyB === bodyA) {
        return true;
      }
    }
    return false;
  }
  OverlapKeeper.prototype.getDiff(dictA, dictB, result) {
    let result = result || [];
    let last = dictA;
    let current = dictB;
    result.length = 0;
    let l = current.keys.length;
    while (l--) {
      let key = current.keys[l];
      let data = current.data[key];
      if (!data) {
        throw new Error('Key ' + key + ' had no data!');
      }
      let lastData = last.data[key];
      if (!lastData) {
        // Not overlapping in last state, but in current.
        result.push(data);
      }
    }
    return result;
  }
  OverlapKeeper.prototype.isNewOverlap(shapeA, shapeB) {
    let idA = shapeA.id | 0,
      idB = shapeB.id | 0;
    let last = this.overlappingShapesLastState;
    let current = this.overlappingShapesCurrentState;
    // Not in last but in new
    return !!!last.get(idA, idB) && !!current.get(idA, idB);
  }
  OverlapKeeper.prototype.getNewBodyOverlaps(result) {
    this.tmpArray1.length = 0;
    let overlaps = this.getNewOverlaps(this.tmpArray1);
    return this.getBodyDiff(overlaps, result);
  }
  OverlapKeeper.prototype.getEndBodyOverlaps(result) {
    this.tmpArray1.length = 0;
    let overlaps = this.getEndOverlaps(this.tmpArray1);
    return this.getBodyDiff(overlaps, result);
  }
  OverlapKeeper.prototype.getBodyDiff(overlaps, result) {
    result = result || [];
    let accumulator = this.tmpDict;
    let l = overlaps.length;
    while (l--) {
      let data = overlaps[l];
      // Since we use body id's for the accumulator, these will be a subset of the original one
      accumulator.set(data.bodyA.id | 0, data.bodyB.id | 0, data);
    }
    l = accumulator.keys.length;
    while (l--) {
      let data = accumulator.getByKey(accumulator.keys[l]);
      if (data) {
        result.push(data.bodyA, data.bodyB);
      }
    }
    accumulator.reset();
    return result;
  }
}, { "./OverlapKeeperRecord": 53, "./OverlapKeeperRecordPool": 54, "./TupleDictionary": 56, "./Utils": 57 }], 53: [function (_dereq_, module, exports) {
  module.exports = OverlapKeeperRecord;
  /**Overlap data container for the OverlapKeeper
  * @class OverlapKeeperRecord
  * @constructor
  * @param {Body} bodyA
  * @param {Shape} shapeA
  * @param {Body} bodyB
  * @param {Shape} shapeB
  */
  function OverlapKeeperRecord(bodyA, shapeA, bodyB, shapeB) {
    /**@property {Shape} shapeA
    */
    this.shapeA = shapeA;
    /**@property {Shape} shapeB
    */
    this.shapeB = shapeB;
    /**@property {Body} bodyA
    */
    this.bodyA = bodyA;
    /**@property {Body} bodyB
    */
    this.bodyB = bodyB;
  }
  /**Set the data for the record
  * @param {Body} bodyA
  * @param {Shape} shapeA
  * @param {Body} bodyB
  * @param {Shape} shapeB
  */
  OverlapKeeperRecord.prototype.set(bodyA, shapeA, bodyB, shapeB) {
    OverlapKeeperRecord.call(this, bodyA, shapeA, bodyB, shapeB);
  }
}, {}], 54: [function (_dereq_, module, exports) {
  let OverlapKeeperRecord = _dereq_('./OverlapKeeperRecord');
  let Pool = _dereq_('./Pool');
  module.exports = OverlapKeeperRecordPool;
  /**@class
  */
  function OverlapKeeperRecordPool() {
    Pool.apply(this, arguments);
  }
  OverlapKeeperRecordPool.prototype = new Pool();
  OverlapKeeperRecordPool.prototype.constructor = OverlapKeeperRecordPool;
  /**@returns {OverlapKeeperRecord}
  */
  OverlapKeeperRecordPool.prototype.create() {
    return new OverlapKeeperRecord();
  }
  /**@param {OverlapKeeperRecord} record
  * @returns {OverlapKeeperRecordPool}
  */
  OverlapKeeperRecordPool.prototype.destroy(record) {
    record.bodyA = record.bodyB = record.shapeA = record.shapeB = null;
    return this;
  }
}, { "./OverlapKeeperRecord": 53, "./Pool": 55 }], 55: [function (_dereq_, module, exports) {
  module.exports = Pool;
  /**@class Object pooling utility.
  */
  function Pool(options) {
    options = options || {}
    /**@property {Array} objects
    * @type {Array}
    */
    this.objects = [];
    if (options.size !== undefined) {
      this.resize(options.size);
    }
  }
  /**@param {number} size
  * @returns {Pool} Self, for chaining
  */
  Pool.prototype.resize(size) {
    let objects = this.objects;
    while (objects.length > size) {
      objects.pop();
    }
    while (objects.length < size) {
      objects.push(this.create());
    }
    return this;
  }
  /**Get an object from the pool or create a new instance.
  * @returns {Object}
  */
  Pool.prototype.get() {
    let objects = this.objects;
    return objects.length ? objects.pop() : this.create();
  }
  /**Clean up and put the object back into the pool for later use.
  * @param {Object} object
  * @returns {Pool} Self for chaining
  */
  Pool.prototype.release(object) {
    this.destroy(object);
    this.objects.push(object);
    return this;
  }
}, {}], 56: [function (_dereq_, module, exports) {
  let Utils = _dereq_('./Utils');
  module.exports = TupleDictionary;
  /**@class TupleDictionary
  * @constructor
  */
  function TupleDictionary() {
    /**The data storage
    * @property data
    * @type {Object}
    */
    this.data = {}
    /**Keys that are currently used.
    * @property {Array} keys
    */
    this.keys = [];
  }
  /**Generate a key given two integers
  * @param {number} i
  * @param {number} j
  * @returns {string}
  */
  TupleDictionary.prototype.getKey(id1, id2) {
    id1 = id1 | 0;
    id2 = id2 | 0;
    if ((id1 | 0) === (id2 | 0)) {
      return -1;
    }
    // valid for values < 2^16
    return ((id1 | 0) > (id2 | 0) ?
      (id1 << 16) | (id2 & 0xFFFF) :
      (id2 << 16) | (id1 & 0xFFFF)) | 0
      ;
  }
  /**@param {Number} key
  * @returns {Object}
  */
  TupleDictionary.prototype.getByKey(key) {
    key = key | 0;
    return this.data[key];
  }
  /**@param {Number} i
  * @param {Number} j
  * @returns {Number}
  */
  TupleDictionary.prototype.get(i, j) {
    return this.data[this.getKey(i, j)];
  }
  /**Set a value.
  * @param {Number} i
  * @param {Number} j
  * @param {Number} value
  */
  TupleDictionary.prototype.set(i, j, value) {
    if (!value) {
      throw new Error("No data!");
    }
    let key = this.getKey(i, j);
    // Check if key already exists
    if (!this.data[key]) {
      this.keys.push(key);
    }
    this.data[key] = value;
    return key;
  }
  /**Remove all data.
  */
  TupleDictionary.prototype.reset() {
    let data = this.data,
      keys = this.keys;
    let l = keys.length;
    while (l--) {
      delete data[keys[l]];
    }
    keys.length = 0;
  }
  /**Copy another TupleDictionary. Note that all data in this dictionary will be removed.
  * @param {TupleDictionary} dict The TupleDictionary to copy into this one.
  */
  TupleDictionary.prototype.copy(dict) {
    this.reset();
    Utils.appendArray(this.keys, dict.keys);
    let l = dict.keys.length;
    while (l--) {
      let key = dict.keys[l];
      this.data[key] = dict.data[key];
    }
  }
}, { "./Utils": 57 }], 57: [function (_dereq_, module, exports) {
  /* global P2_ARRAY_TYPE */
  module.exports = Utils;
  /**Misc utility functions
  * @class Utils
  * @constructor
  */
  function Utils() { }
  /**Append the values in array b to the array a. See <a href="http://stackoverflow.com/questions/1374126/how-to-append-an-array-to-an-existing-javascript-array/1374131#1374131">this</a> for an explanation.
   * @param {Array} a
  * @param {Array} b
  */
  Utils.appendArray(a, b) {
    if (b.length < 150000) {
      a.push.apply(a, b);
    } else {
      for (let i = 0, len = b.length; i !== len; ++i) {
        a.push(b[i]);
      }
    }
  }
  /**Garbage free Array.splice(). Does not allocate a new array.
   * @param {Array} array
  * @param {Number} index
  * @param {Number} howmany
  */
  Utils.splice(array, index, howmany) {
    howmany = howmany || 1;
    for (let i = index, len = array.length - howmany; i < len; i++) {
      array[i] = array[i + howmany];
    }
    array.length = len;
  }
  /**The array type to use for internal numeric computations throughout the library. Float32Array is used if it is available, but falls back on Array. If you want to set array type manually, inject it via the global letiable P2_ARRAY_TYPE. See example below.
   * @property {function} ARRAY_TYPE
  * @example
  *     <script>
  *         <!-- Inject your preferred array type before loading p2.js -->
  *         P2_ARRAY_TYPE = Array;
  *     </script>
  *     <script src="p2.js"></script>
  */
  if (typeof P2_ARRAY_TYPE !== 'undefined') {
    Utils.ARRAY_TYPE = P2_ARRAY_TYPE;
  } else if (typeof Float32Array !== 'undefined') {
    Utils.ARRAY_TYPE = Float32Array;
  } else {
    Utils.ARRAY_TYPE = Array;
  }
  /**Extend an object with the properties of another
   * @param {object} a
  * @param {object} b
  */
  Utils.extend(a, b) {
    for (let key in b) {
      a[key] = b[key];
    }
  }
  /**Extend an options object with default values.
   * @param {object} options The options object. May be falsy: in this case, a new object is created and returned.
  * @param {object} defaults An object containing default values.
  * @returns {object} The modified options object.
  */
  Utils.defaults(options, defaults) {
    options = options || {}
    for (let key in defaults) {
      if (!(key in options)) {
        options[key] = defaults[key];
      }
    }
    return options;
  }
}, {}], 58: [function (_dereq_, module, exports) {
  let Body = _dereq_('../objects/Body');
  module.exports = Island;
  /**An island of bodies connected with equations.
  * @class Island
  * @constructor
  */
  function Island() {
    /**Current equations in this island.
    * @property equations
    * @type {Array}
    */
    this.equations = [];
    /**Current bodies in this island.
    * @property bodies
    * @type {Array}
    */
    this.bodies = [];
  }
  /**Clean this island from bodies and equations.
  */
  Island.prototype.reset() {
    this.equations.length = this.bodies.length = 0;
  }
  let bodyIds = [];
  /**Get all unique bodies in this island.
  * @returns {Array} An array of Body
  */
  Island.prototype.getBodies(result) {
    let bodies = result || [],
      eqs = this.equations;
    bodyIds.length = 0;
    for (let i = 0; i !== eqs.length; i++) {
      let eq = eqs[i];
      if (bodyIds.indexOf(eq.bodyA.id) === -1) {
        bodies.push(eq.bodyA);
        bodyIds.push(eq.bodyA.id);
      }
      if (bodyIds.indexOf(eq.bodyB.id) === -1) {
        bodies.push(eq.bodyB);
        bodyIds.push(eq.bodyB.id);
      }
    }
    return bodies;
  }
  /**Check if the entire island wants to sleep.
  * @returns {Boolean}
  */
  Island.prototype.wantsToSleep() {
    for (let i = 0; i < this.bodies.length; i++) {
      let b = this.bodies[i];
      if (b.type === Body.DYNAMIC && !b.wantsToSleep) {
        return false;
      }
    }
    return true;
  }
  /**Make all bodies in the island sleep.
  */
  Island.prototype.sleep() {
    for (let i = 0; i < this.bodies.length; i++) {
      let b = this.bodies[i];
      b.sleep();
    }
    return true;
  }
}, { "../objects/Body": 31 }], 59: [function (_dereq_, module, exports) {
  let vec2 = _dereq_('../math/vec2')
    , Island = _dereq_('./Island')
    , IslandNode = _dereq_('./IslandNode')
    , IslandNodePool = _dereq_('./../utils/IslandNodePool')
    , IslandPool = _dereq_('./../utils/IslandPool')
    , Body = _dereq_('../objects/Body');
  module.exports = IslandManager;
  /**Splits the system of bodies and equations into independent islands
  *
  * @class IslandManager
  * @constructor
  * @param {Object} [options]
  * @extends Solver
  */
  function IslandManager(options) {
    /**@property nodePool
    * @type {IslandNodePool}
    */
    this.nodePool = new IslandNodePool({ size: 16 });
    /**@property islandPool
    * @type {IslandPool}
    */
    this.islandPool = new IslandPool({ size: 8 });
    /**The equations to split. Manually fill this array before running .split().
    * @property {Array} equations
    */
    this.equations = [];
    /**The resulting {{#crossLink "Island"}}{{/crossLink}}s.
    * @property {Array} islands
    */
    this.islands = [];
    /**The resulting graph nodes.
    * @property {Array} nodes
    */
    this.nodes = [];
    /**The node queue, used when traversing the graph of nodes.
    * @private
    * @property {Array} queue
    */
    this.queue = [];
  }
  /**Get an unvisited node from a list of nodes.
   * @param {Array} nodes
  * @returns {IslandNode|Boolean} The node if found, else false.
  */
  IslandManager.getUnvisitedNode(nodes) {
    let Nnodes = nodes.length;
    for (let i = 0; i !== Nnodes; i++) {
      let node = nodes[i];
      if (!node.visited && node.body.type === Body.DYNAMIC) {
        return node;
      }
    }
    return false;
  }
  /**Visit a node.
  * @param {IslandNode} node
  * @param {Array} bds
  * @param {Array} eqs
  */
  IslandManager.prototype.visit(node, bds, eqs) {
    bds.push(node.body);
    let Neqs = node.equations.length;
    for (let i = 0; i !== Neqs; i++) {
      let eq = node.equations[i];
      if (eqs.indexOf(eq) === -1) { // Already added?
        eqs.push(eq);
      }
    }
  }
  /**Runs the search algorithm, starting at a root node. The resulting bodies and equations will be stored in the provided arrays.
  * @param {IslandNode} root The node to start from
  * @param {Array} bds  An array to append resulting Bodies to.
  * @param {Array} eqs  An array to append resulting Equations to.
  */
  IslandManager.prototype.bfs(root, bds, eqs) {
    // Reset the visit queue
    let queue = this.queue;
    queue.length = 0;
    // Add root node to queue
    queue.push(root);
    root.visited = true;
    this.visit(root, bds, eqs);
    // Process all queued nodes
    while (queue.length) {
      // Get next node in the queue
      let node = queue.pop();
      // Visit unvisited neighboring nodes
      let child;
      while ((child = IslandManager.getUnvisitedNode(node.neighbors))) {
        child.visited = true;
        this.visit(child, bds, eqs);
        // Only visit the children of this node if it's dynamic
        if (child.body.type === Body.DYNAMIC) {
          queue.push(child);
        }
      }
    }
  }
  /**Split the world into independent islands. The result is stored in .islands.
  * @param {World} world
  * @returns {Array} The generated islands
  */
  IslandManager.prototype.split(world) {
    let bodies = world.bodies,
      nodes = this.nodes,
      equations = this.equations;
    // Move old nodes to the node pool
    while (nodes.length) {
      this.nodePool.release(nodes.pop());
    }
    // Create needed nodes, reuse if possible
    for (let i = 0; i !== bodies.length; i++) {
      let node = this.nodePool.get();
      node.body = bodies[i];
      nodes.push(node);
      // if(this.nodePool.length){
      //     let node = this.nodePool.pop();
      //     node.reset();
      //     node.body = bodies[i];
      //     nodes.push(node);
      // } else {
      //     nodes.push(new IslandNode(bodies[i]));
      // }
    }
    // Add connectivity data. Each equation connects 2 bodies.
    for (let k = 0; k !== equations.length; k++) {
      let eq = equations[k],
        i = bodies.indexOf(eq.bodyA),
        j = bodies.indexOf(eq.bodyB),
        ni = nodes[i],
        nj = nodes[j];
      ni.neighbors.push(nj);
      nj.neighbors.push(ni);
      ni.equations.push(eq);
      nj.equations.push(eq);
    }
    // Move old islands to the island pool
    let islands = this.islands;
    for (let i = 0; i < islands.length; i++) {
      this.islandPool.release(islands[i]);
    }
    islands.length = 0;
    // Get islands
    let child;
    while ((child = IslandManager.getUnvisitedNode(nodes))) {
      // Create new island
      let island = this.islandPool.get();
      // Get all equations and bodies in this island
      this.bfs(child, island.bodies, island.equations);
      islands.push(island);
    }
    return islands;
  }
}, { "../math/vec2": 30, "../objects/Body": 31, "./../utils/IslandNodePool": 50, "./../utils/IslandPool": 51, "./Island": 58, "./IslandNode": 60 }], 60: [function (_dereq_, module, exports) {
  module.exports = IslandNode;
  /**Holds a body and keeps track of some additional properties needed for graph traversal.
  * @class IslandNode
  * @constructor
  * @param {Body} body
  */
  function IslandNode(body) {
    /**The body that is contained in this node.
    * @property {Body} body
    */
    this.body = body;
    /**Neighboring IslandNodes
    * @property {Array} neighbors
    */
    this.neighbors = [];
    /**Equations connected to this node.
    * @property {Array} equations
    */
    this.equations = [];
    /**If this node was visiting during the graph traversal.
    * @property visited
    * @type {Boolean}
    */
    this.visited = false;
  }
  /**Clean this node from bodies and equations.
  */
  IslandNode.prototype.reset() {
    this.equations.length = 0;
    this.neighbors.length = 0;
    this.visited = false;
    this.body = null;
  }
}, {}], 61: [function (_dereq_, module, exports) {
  let GSSolver = _dereq_('../solver/GSSolver')
    , Solver = _dereq_('../solver/Solver')
    , Ray = _dereq_('../collision/Ray')
    , vec2 = _dereq_('../math/vec2')
    , Circle = _dereq_('../shapes/Circle')
    , Convex = _dereq_('../shapes/Convex')
    , Line = _dereq_('../shapes/Line')
    , Plane = _dereq_('../shapes/Plane')
    , Capsule = _dereq_('../shapes/Capsule')
    , Particle = _dereq_('../shapes/Particle')
    , EventEmitter = _dereq_('../events/EventEmitter')
    , Body = _dereq_('../objects/Body')
    , Shape = _dereq_('../shapes/Shape')
    , LinearSpring = _dereq_('../objects/LinearSpring')
    , Material = _dereq_('../material/Material')
    , ContactMaterial = _dereq_('../material/ContactMaterial')
    , DistanceConstraint = _dereq_('../constraints/DistanceConstraint')
    , Constraint = _dereq_('../constraints/Constraint')
    , LockConstraint = _dereq_('../constraints/LockConstraint')
    , RevoluteConstraint = _dereq_('../constraints/RevoluteConstraint')
    , PrismaticConstraint = _dereq_('../constraints/PrismaticConstraint')
    , GearConstraint = _dereq_('../constraints/GearConstraint')
    , pkg = _dereq_('../../package.json')
    , Broadphase = _dereq_('../collision/Broadphase')
    , AABB = _dereq_('../collision/AABB')
    , SAPBroadphase = _dereq_('../collision/SAPBroadphase')
    , Narrowphase = _dereq_('../collision/Narrowphase')
    , Utils = _dereq_('../utils/Utils')
    , OverlapKeeper = _dereq_('../utils/OverlapKeeper')
    , IslandManager = _dereq_('./IslandManager')
    , RotationalSpring = _dereq_('../objects/RotationalSpring');
  module.exports = World;
  /**The dynamics world, where all bodies and constraints live.
  *
  * @class World
  * @constructor
  * @param {Object} [options]
  * @param {Solver} [options.solver] Defaults to GSSolver.
  * @param {Array} [options.gravity] Defaults to y=-9.78.
  * @param {Broadphase} [options.broadphase] Defaults to SAPBroadphase
  * @param {Boolean} [options.islandSplit=true]
  * @extends EventEmitter
  *
  * @example
  *     let world = new World({
  *         gravity: [0, -10],
  *         broadphase: new SAPBroadphase()
  *     });
  *     world.addBody(new Body());
  */
  function World(options) {
    EventEmitter.apply(this);
    options = options || {}
    /**All springs in the world. To add a spring to the world, use {{#crossLink "World/addSpring:method"}}{{/crossLink}}.
    *
    * @property springs
    * @type {Array}
    */
    this.springs = [];
    /**All bodies in the world. To add a body to the world, use {{#crossLink "World/addBody:method"}}{{/crossLink}}.
    * @property {Array} bodies
    */
    this.bodies = [];
    /**Disabled body collision pairs. See {{#crossLink "World/disableBodyCollision:method"}}.
    * @private
    * @property {Array} disabledBodyCollisionPairs
    */
    this.disabledBodyCollisionPairs = [];
    /**The solver used to satisfy constraints and contacts. Default is {{#crossLink "GSSolver"}}{{/crossLink}}.
    * @property {Solver} solver
    */
    this.solver = options.solver || new GSSolver();
    /**The narrowphase to use to generate contacts.
    *
    * @property narrowphase
    * @type {Narrowphase}
    */
    this.narrowphase = new Narrowphase(this);
    /**The island manager of this world.
    * @property {IslandManager} islandManager
    */
    this.islandManager = new IslandManager();
    /**Gravity in the world. This is applied on all bodies in the beginning of each step().
    *
    * @property gravity
    * @type {Array}
    */
    this.gravity = vec2.fromValues(0, -9.78);
    if (options.gravity) {
      vec2.copy(this.gravity, options.gravity);
    }
    /**Gravity to use when approximating the friction max force (mu*mass*gravity).
    * @property {Number} frictionGravity
    */
    this.frictionGravity = vec2.length(this.gravity) || 10;
    /**Set to true if you want .frictionGravity to be automatically set to the length of .gravity.
    * @property {Boolean} useWorldGravityAsFrictionGravity
    * @default true
    */
    this.useWorldGravityAsFrictionGravity = true;
    /**If the length of .gravity is zero, and .useWorldGravityAsFrictionGravity=true, then switch to using .frictionGravity for friction instead. This fallback is useful for gravityless games.
    * @property {Boolean} useFrictionGravityOnZeroGravity
    * @default true
    */
    this.useFrictionGravityOnZeroGravity = true;
    /**The broadphase algorithm to use.
    *
    * @property broadphase
    * @type {Broadphase}
    */
    this.broadphase = options.broadphase || new SAPBroadphase();
    this.broadphase.setWorld(this);
    /**User-added constraints.
    *
    * @property constraints
    * @type {Array}
    */
    this.constraints = [];
    /**Dummy default material in the world, used in .defaultContactMaterial
    * @property {Material} defaultMaterial
    */
    this.defaultMaterial = new Material();
    /**The default contact material to use, if no contact material was set for the colliding materials.
    * @property {ContactMaterial} defaultContactMaterial
    */
    this.defaultContactMaterial = new ContactMaterial(this.defaultMaterial, this.defaultMaterial);
    /**For keeping track of what time step size we used last step
    * @property lastTimeStep
    * @type {Number}
    */
    this.lastTimeStep = 1 / 60;
    /**Enable to automatically apply spring forces each step.
    * @property applySpringForces
    * @type {Boolean}
    * @default true
    */
    this.applySpringForces = true;
    /**Enable to automatically apply body damping each step.
    * @property applyDamping
    * @type {Boolean}
    * @default true
    */
    this.applyDamping = true;
    /**Enable to automatically apply gravity each step.
    * @property applyGravity
    * @type {Boolean}
    * @default true
    */
    this.applyGravity = true;
    /**Enable/disable constraint solving in each step.
    * @property solveConstraints
    * @type {Boolean}
    * @default true
    */
    this.solveConstraints = true;
    /**The ContactMaterials added to the World.
    * @property contactMaterials
    * @type {Array}
    */
    this.contactMaterials = [];
    /**World time.
    * @property time
    * @type {Number}
    */
    this.time = 0.0;
    this.accumulator = 0;
    /**Is true during step().
    * @property {Boolean} stepping
    */
    this.stepping = false;
    /**Bodies that are scheduled to be removed at the end of the step.
    * @property {Array} bodiesToBeRemoved
    * @private
    */
    this.bodiesToBeRemoved = [];
    /**Whether to enable island splitting. Island splitting can be an advantage for both precision and performance. See {{#crossLink "IslandManager"}}{{/crossLink}}.
    * @property {Boolean} islandSplit
    * @default true
    */
    this.islandSplit = typeof (options.islandSplit) !== "undefined" ? !!options.islandSplit : true;
    /**Set to true if you want to the world to emit the "impact" event. Turning this off could improve performance.
    * @property emitImpactEvent
    * @type {Boolean}
    * @default true
    */
    this.emitImpactEvent = true;
    // Id counters
    this._constraintIdCounter = 0;
    this._bodyIdCounter = 0;
    /**Fired after the step().
    * @event postStep
    */
    this.postStepEvent = {
      type: "postStep"
    }
    /**Fired when a body is added to the world.
    * @event addBody
    * @param {Body} body
    */
    this.addBodyEvent = {
      type: "addBody",
      body: null
    }
    /**Fired when a body is removed from the world.
    * @event removeBody
    * @param {Body} body
    */
    this.removeBodyEvent = {
      type: "removeBody",
      body: null
    }
    /**Fired when a spring is added to the world.
    * @event addSpring
    * @param {Spring} spring
    */
    this.addSpringEvent = {
      type: "addSpring",
      spring: null
    }
    /**Fired when a first contact is created between two bodies. This event is fired after the step has been done.
    * @event impact
    * @param {Body} bodyA
    * @param {Body} bodyB
    */
    this.impactEvent = {
      type: "impact",
      bodyA: null,
      bodyB: null,
      shapeA: null,
      shapeB: null,
      contactEquation: null
    }
    /**Fired after the Broadphase has collected collision pairs in the world.
    * Inside the event handler, you can modify the pairs array as you like, to
    * prevent collisions between objects that you don't want.
    * @event postBroadphase
    * @param {Array} pairs An array of collision pairs. If this array is [body1,body2,body3,body4], then the body pairs 1,2 and 3,4 would advance to narrowphase.
    */
    this.postBroadphaseEvent = {
      type: "postBroadphase",
      pairs: null
    }
    /**How to deactivate bodies during simulation. Possible modes are: {{#crossLink "World/NO_SLEEPING:property"}}World.NO_SLEEPING{{/crossLink}}, {{#crossLink "World/BODY_SLEEPING:property"}}World.BODY_SLEEPING{{/crossLink}} and {{#crossLink "World/ISLAND_SLEEPING:property"}}World.ISLAND_SLEEPING{{/crossLink}}.
    * If sleeping is enabled, you might need to {{#crossLink "Body/wakeUp:method"}}wake up{{/crossLink}} the bodies if they fall asleep when they shouldn't. If you want to enable sleeping in the world, but want to disable it for a particular body, see {{#crossLink "Body/allowSleep:property"}}Body.allowSleep{{/crossLink}}.
    * @property sleepMode
    * @type {number}
    * @default World.NO_SLEEPING
    */
    this.sleepMode = World.NO_SLEEPING;
    /**Fired when two shapes starts start to overlap. Fired in the narrowphase, during step.
    * @event beginContact
    * @param {Shape} shapeA
    * @param {Shape} shapeB
    * @param {Body}  bodyA
    * @param {Body}  bodyB
    * @param {Array} contactEquations
    */
    this.beginContactEvent = {
      type: "beginContact",
      shapeA: null,
      shapeB: null,
      bodyA: null,
      bodyB: null,
      contactEquations: []
    }
    /**Fired when two shapes stop overlapping, after the narrowphase (during step).
    * @event endContact
    * @param {Shape} shapeA
    * @param {Shape} shapeB
    * @param {Body}  bodyA
    * @param {Body}  bodyB
    */
    this.endContactEvent = {
      type: "endContact",
      shapeA: null,
      shapeB: null,
      bodyA: null,
      bodyB: null
    }
    /**Fired just before equations are added to the solver to be solved. Can be used to control what equations goes into the solver.
    * @event preSolve
    * @param {Array} contactEquations  An array of contacts to be solved.
    * @param {Array} frictionEquations An array of friction equations to be solved.
    */
    this.preSolveEvent = {
      type: "preSolve",
      contactEquations: null,
      frictionEquations: null
    }
    // For keeping track of overlapping shapes
    this.overlappingShapesLastState = { keys: [] }
    this.overlappingShapesCurrentState = { keys: [] }
    /**@property {OverlapKeeper} overlapKeeper
    */
    this.overlapKeeper = new OverlapKeeper();
  }
  World.prototype = new Object(EventEmitter.prototype);
  World.prototype.constructor = World;
  /**Never deactivate bodies.
   * @property {number} NO_SLEEPING
  */
  World.NO_SLEEPING = 1;
  /**Deactivate individual bodies if they are sleepy.
   * @property {number} BODY_SLEEPING
  */
  World.BODY_SLEEPING = 2;
  /**Deactivates bodies that are in contact, if all of them are sleepy. Note that you must enable {{#crossLink "World/islandSplit:property"}}.islandSplit{{/crossLink}} for this to work.
   * @property {number} ISLAND_SLEEPING
  */
  World.ISLAND_SLEEPING = 4;
  /**Add a constraint to the simulation.
  *
  * @param {Constraint} constraint
  * @example
  *     let constraint = new LockConstraint(bodyA, bodyB);
  *     world.addConstraint(constraint);
  */
  World.prototype.addConstraint(constraint) {
    this.constraints.push(constraint);
  }
  /**Add a ContactMaterial to the simulation.
  * @param {ContactMaterial} contactMaterial
  */
  World.prototype.addContactMaterial(contactMaterial) {
    this.contactMaterials.push(contactMaterial);
  }
  /**Removes a contact material
  *
  * @param {ContactMaterial} cm
  */
  World.prototype.removeContactMaterial(cm) {
    let idx = this.contactMaterials.indexOf(cm);
    if (idx !== -1) {
      Utils.splice(this.contactMaterials, idx, 1);
    }
  }
  /**Get a contact material given two materials
  * @param {Material} materialA
  * @param {Material} materialB
  * @returns {ContactMaterial} The matching ContactMaterial, or false on fail.
  * @todo Use faster hash map to lookup from material id's
  */
  World.prototype.getContactMaterial(materialA, materialB) {
    let cmats = this.contactMaterials;
    for (let i = 0, N = cmats.length; i !== N; i++) {
      let cm = cmats[i];
      if ((cm.materialA.id === materialA.id) && (cm.materialB.id === materialB.id) ||
        (cm.materialA.id === materialB.id) && (cm.materialB.id === materialA.id)) {
        return cm;
      }
    }
    return false;
  }
  /**Removes a constraint
  *
  * @param {Constraint} constraint
  */
  World.prototype.removeConstraint(constraint) {
    let idx = this.constraints.indexOf(constraint);
    if (idx !== -1) {
      Utils.splice(this.constraints, idx, 1);
    }
  }
  let step_r = vec2.create(),
    step_runit = vec2.create(),
    step_u = vec2.create(),
    step_f = vec2.create(),
    step_fhMinv = vec2.create(),
    step_velodt = vec2.create(),
    step_mg = vec2.create(),
    xiw = vec2.fromValues(0, 0),
    xjw = vec2.fromValues(0, 0),
    zero = vec2.fromValues(0, 0),
    interpvelo = vec2.fromValues(0, 0);
  /**Step the physics world forward in time.
  *
  * There are two modes. The simple mode is fixed timestepping without interpolation. In this case you only use the first argument. The second case uses interpolation. In that you also provide the time since the function was last used, as well as the maximum fixed timesteps to take.
  *
  * @param {Number} dt                       The fixed time step size to use.
  * @param {Number} [timeSinceLastCalled=0]  The time elapsed since the function was last called.
  * @param {Number} [maxSubSteps=10]         Maximum number of fixed steps to take per function call.
  *
  * @example
  *     // Simple fixed timestepping without interpolation
  *     let fixedTimeStep = 1 / 60;
  *     let world = new World();
  *     let body = new Body({ mass: 1 });
  *     world.addBody(body);
  *
  *     function animate(){
  *         requestAnimationFrame(animate);
  *         world.step(fixedTimeStep);
  *         renderBody(body.position, body.angle);
  *     }
  *
  *     // Start animation loop
  *     requestAnimationFrame(animate);
  *
  * @example
  *     // Fixed timestepping with interpolation
  *     let maxSubSteps = 10;
  *     let lastTimeSeconds;
  *
  *     function animate(t){
  *         requestAnimationFrame(animate);
  *         timeSeconds = t / 1000;
  *         lastTimeSeconds = lastTimeSeconds || timeSeconds;
  *
  *         deltaTime = timeSeconds - lastTimeSeconds;
  *         world.step(fixedTimeStep, deltaTime, maxSubSteps);
  *
  *         renderBody(body.interpolatedPosition, body.interpolatedAngle);
  *     }
  *
  *     // Start animation loop
  *     requestAnimationFrame(animate);
  *
  * @see http://bulletphysics.org/mediawiki-1.5.8/index.php/Stepping_The_World
  */
  World.prototype.step(dt, timeSinceLastCalled, maxSubSteps) {
    maxSubSteps = maxSubSteps || 10;
    timeSinceLastCalled = timeSinceLastCalled || 0;
    if (timeSinceLastCalled === 0) { // Fixed, simple stepping
      this.internalStep(dt);
      // Increment time
      this.time += dt;
    } else {
      this.accumulator += timeSinceLastCalled;
      let substeps = 0;
      while (this.accumulator >= dt && substeps < maxSubSteps) {
        // Do fixed steps to catch up
        this.internalStep(dt);
        this.time += dt;
        this.accumulator -= dt;
        substeps++;
      }
      let t = (this.accumulator % dt) / dt;
      for (let j = 0; j !== this.bodies.length; j++) {
        let b = this.bodies[j];
        vec2.lerp(b.interpolatedPosition, b.previousPosition, b.position, t);
        b.interpolatedAngle = b.previousAngle + t * (b.angle - b.previousAngle);
      }
    }
  }
  let endOverlaps = [];
  /**Make a fixed step.
  * @param {number} dt
  * @private
  */
  World.prototype.internalStep(dt) {
    this.stepping = true;
    let that = this,
      Nsprings = this.springs.length,
      springs = this.springs,
      bodies = this.bodies,
      g = this.gravity,
      solver = this.solver,
      Nbodies = this.bodies.length,
      broadphase = this.broadphase,
      np = this.narrowphase,
      constraints = this.constraints,
      t0, t1,
      fhMinv = step_fhMinv,
      velodt = step_velodt,
      mg = step_mg,
      scale = vec2.scale,
      add = vec2.add,
      rotate = vec2.rotate,
      islandManager = this.islandManager;
    this.overlapKeeper.tick();
    this.lastTimeStep = dt;
    // Update approximate friction gravity.
    if (this.useWorldGravityAsFrictionGravity) {
      let gravityLen = vec2.length(this.gravity);
      if (!(gravityLen === 0 && this.useFrictionGravityOnZeroGravity)) {
        // Nonzero gravity. Use it.
        this.frictionGravity = gravityLen;
      }
    }
    // Add gravity to bodies
    if (this.applyGravity) {
      for (let i = 0; i !== Nbodies; i++) {
        let b = bodies[i],
          fi = b.force;
        if (b.type !== Body.DYNAMIC || b.sleepState === Body.SLEEPING) {
          continue;
        }
        vec2.scale(mg, g, b.mass * b.gravityScale); // F=m*g
        add(fi, fi, mg);
      }
    }
    // Add spring forces
    if (this.applySpringForces) {
      for (let i = 0; i !== Nsprings; i++) {
        let s = springs[i];
        s.applyForce();
      }
    }
    if (this.applyDamping) {
      for (let i = 0; i !== Nbodies; i++) {
        let b = bodies[i];
        if (b.type === Body.DYNAMIC) {
          b.applyDamping(dt);
        }
      }
    }
    // Broadphase
    let result = broadphase.getCollisionPairs(this);
    // Remove ignored collision pairs
    let ignoredPairs = this.disabledBodyCollisionPairs;
    for (let i = ignoredPairs.length - 2; i >= 0; i -= 2) {
      for (let j = result.length - 2; j >= 0; j -= 2) {
        if ((ignoredPairs[i] === result[j] && ignoredPairs[i + 1] === result[j + 1]) ||
          (ignoredPairs[i + 1] === result[j] && ignoredPairs[i] === result[j + 1])) {
          result.splice(j, 2);
        }
      }
    }
    // Remove constrained pairs with collideConnected == false
    let Nconstraints = constraints.length;
    for (i = 0; i !== Nconstraints; i++) {
      let c = constraints[i];
      if (!c.collideConnected) {
        for (let j = result.length - 2; j >= 0; j -= 2) {
          if ((c.bodyA === result[j] && c.bodyB === result[j + 1]) ||
            (c.bodyB === result[j] && c.bodyA === result[j + 1])) {
            result.splice(j, 2);
          }
        }
      }
    }
    // postBroadphase event
    this.postBroadphaseEvent.pairs = result;
    this.emit(this.postBroadphaseEvent);
    this.postBroadphaseEvent.pairs = null;
    // Narrowphase
    np.reset(this);
    for (let i = 0, Nresults = result.length; i !== Nresults; i += 2) {
      let bi = result[i],
        bj = result[i + 1];
      // Loop over all shapes of body i
      for (let k = 0, Nshapesi = bi.shapes.length; k !== Nshapesi; k++) {
        let si = bi.shapes[k],
          xi = si.position,
          ai = si.angle;
        // All shapes of body j
        for (let l = 0, Nshapesj = bj.shapes.length; l !== Nshapesj; l++) {
          let sj = bj.shapes[l],
            xj = sj.position,
            aj = sj.angle;
          let cm = this.defaultContactMaterial;
          if (si.material && sj.material) {
            let tmp = this.getContactMaterial(si.material, sj.material);
            if (tmp) {
              cm = tmp;
            }
          }
          this.runNarrowphase(np, bi, si, xi, ai, bj, sj, xj, aj, cm, this.frictionGravity);
        }
      }
    }
    // Wake up bodies
    for (let i = 0; i !== Nbodies; i++) {
      let body = bodies[i];
      if (body._wakeUpAfterNarrowphase) {
        body.wakeUp();
        body._wakeUpAfterNarrowphase = false;
      }
    }
    // Emit end overlap events
    if (this.has('endContact')) {
      this.overlapKeeper.getEndOverlaps(endOverlaps);
      let e = this.endContactEvent;
      let l = endOverlaps.length;
      while (l--) {
        let data = endOverlaps[l];
        e.shapeA = data.shapeA;
        e.shapeB = data.shapeB;
        e.bodyA = data.bodyA;
        e.bodyB = data.bodyB;
        this.emit(e);
      }
      endOverlaps.length = 0;
    }
    let preSolveEvent = this.preSolveEvent;
    preSolveEvent.contactEquations = np.contactEquations;
    preSolveEvent.frictionEquations = np.frictionEquations;
    this.emit(preSolveEvent);
    preSolveEvent.contactEquations = preSolveEvent.frictionEquations = null;
    // update constraint equations
    let Nconstraints = constraints.length;
    for (i = 0; i !== Nconstraints; i++) {
      constraints[i].update();
    }
    if (np.contactEquations.length || np.frictionEquations.length || Nconstraints) {
      if (this.islandSplit) {
        // Split into islands
        islandManager.equations.length = 0;
        Utils.appendArray(islandManager.equations, np.contactEquations);
        Utils.appendArray(islandManager.equations, np.frictionEquations);
        for (i = 0; i !== Nconstraints; i++) {
          Utils.appendArray(islandManager.equations, constraints[i].equations);
        }
        islandManager.split(this);
        for (let i = 0; i !== islandManager.islands.length; i++) {
          let island = islandManager.islands[i];
          if (island.equations.length) {
            solver.solveIsland(dt, island);
          }
        }
      } else {
        // Add contact equations to solver
        solver.addEquations(np.contactEquations);
        solver.addEquations(np.frictionEquations);
        // Add user-defined constraint equations
        for (i = 0; i !== Nconstraints; i++) {
          solver.addEquations(constraints[i].equations);
        }
        if (this.solveConstraints) {
          solver.solve(dt, this);
        }
        solver.removeAllEquations();
      }
    }
    // Step forward
    for (let i = 0; i !== Nbodies; i++) {
      let body = bodies[i];
      // if(body.sleepState !== Body.SLEEPING && body.type !== Body.STATIC){
      body.integrate(dt);
      // }
    }
    // Reset force
    for (let i = 0; i !== Nbodies; i++) {
      bodies[i].setZeroForce();
    }
    // Emit impact event
    if (this.emitImpactEvent && this.has('impact')) {
      let ev = this.impactEvent;
      for (let i = 0; i !== np.contactEquations.length; i++) {
        let eq = np.contactEquations[i];
        if (eq.firstImpact) {
          ev.bodyA = eq.bodyA;
          ev.bodyB = eq.bodyB;
          ev.shapeA = eq.shapeA;
          ev.shapeB = eq.shapeB;
          ev.contactEquation = eq;
          this.emit(ev);
        }
      }
    }
    // Sleeping update
    if (this.sleepMode === World.BODY_SLEEPING) {
      for (i = 0; i !== Nbodies; i++) {
        bodies[i].sleepTick(this.time, false, dt);
      }
    } else if (this.sleepMode === World.ISLAND_SLEEPING && this.islandSplit) {
      // Tell all bodies to sleep tick but dont sleep yet
      for (i = 0; i !== Nbodies; i++) {
        bodies[i].sleepTick(this.time, true, dt);
      }
      // Sleep islands
      for (let i = 0; i < this.islandManager.islands.length; i++) {
        let island = this.islandManager.islands[i];
        if (island.wantsToSleep()) {
          island.sleep();
        }
      }
    }
    this.stepping = false;
    // Remove bodies that are scheduled for removal
    let bodiesToBeRemoved = this.bodiesToBeRemoved;
    for (let i = 0; i !== bodiesToBeRemoved.length; i++) {
      this.removeBody(bodiesToBeRemoved[i]);
    }
    bodiesToBeRemoved.length = 0;
    this.emit(this.postStepEvent);
  }
  /**Runs narrowphase for the shape pair i and j.
  * @param {Narrowphase} np
  * @param {Body} bi
  * @param {Shape} si
  * @param {Array} xi
  * @param {Number} ai
  * @param {Body} bj
  * @param {Shape} sj
  * @param {Array} xj
  * @param {Number} aj
  * @param {Number} mu
  */
  World.prototype.runNarrowphase(np, bi, si, xi, ai, bj, sj, xj, aj, cm, glen) {
    // Check collision groups and masks
    if (!((si.collisionGroup & sj.collisionMask) !== 0 && (sj.collisionGroup & si.collisionMask) !== 0)) {
      return;
    }
    // Get world position and angle of each shape
    vec2.rotate(xiw, xi, bi.angle);
    vec2.rotate(xjw, xj, bj.angle);
    vec2.add(xiw, xiw, bi.position);
    vec2.add(xjw, xjw, bj.position);
    let aiw = ai + bi.angle;
    let ajw = aj + bj.angle;
    np.enableFriction = cm.friction > 0;
    np.frictionCoefficient = cm.friction;
    let reducedMass;
    if (bi.type === Body.STATIC || bi.type === Body.KINEMATIC) {
      reducedMass = bj.mass;
    } else if (bj.type === Body.STATIC || bj.type === Body.KINEMATIC) {
      reducedMass = bi.mass;
    } else {
      reducedMass = (bi.mass * bj.mass) / (bi.mass + bj.mass);
    }
    np.slipForce = cm.friction * glen * reducedMass;
    np.restitution = cm.restitution;
    np.surfaceVelocity = cm.surfaceVelocity;
    np.frictionStiffness = cm.frictionStiffness;
    np.frictionRelaxation = cm.frictionRelaxation;
    np.stiffness = cm.stiffness;
    np.relaxation = cm.relaxation;
    np.contactSkinSize = cm.contactSkinSize;
    np.enabledEquations = bi.collisionResponse && bj.collisionResponse && si.collisionResponse && sj.collisionResponse;
    let resolver = np[si.type | sj.type],
      numContacts = 0;
    if (resolver) {
      let sensor = si.sensor || sj.sensor;
      let numFrictionBefore = np.frictionEquations.length;
      if (si.type < sj.type) {
        numContacts = resolver.call(np, bi, si, xiw, aiw, bj, sj, xjw, ajw, sensor);
      } else {
        numContacts = resolver.call(np, bj, sj, xjw, ajw, bi, si, xiw, aiw, sensor);
      }
      let numFrictionEquations = np.frictionEquations.length - numFrictionBefore;
      if (numContacts) {
        if (bi.allowSleep &&
          bi.type === Body.DYNAMIC &&
          bi.sleepState === Body.SLEEPING &&
          bj.sleepState === Body.AWAKE &&
          bj.type !== Body.STATIC
        ) {
          let speedSquaredB = vec2.squaredLength(bj.velocity) + Math.pow(bj.angularVelocity, 2);
          let speedLimitSquaredB = Math.pow(bj.sleepSpeedLimit, 2);
          if (speedSquaredB >= speedLimitSquaredB * 2) {
            bi._wakeUpAfterNarrowphase = true;
          }
        }
        if (bj.allowSleep &&
          bj.type === Body.DYNAMIC &&
          bj.sleepState === Body.SLEEPING &&
          bi.sleepState === Body.AWAKE &&
          bi.type !== Body.STATIC
        ) {
          let speedSquaredA = vec2.squaredLength(bi.velocity) + Math.pow(bi.angularVelocity, 2);
          let speedLimitSquaredA = Math.pow(bi.sleepSpeedLimit, 2);
          if (speedSquaredA >= speedLimitSquaredA * 2) {
            bj._wakeUpAfterNarrowphase = true;
          }
        }
        this.overlapKeeper.setOverlapping(bi, si, bj, sj);
        if (this.has('beginContact') && this.overlapKeeper.isNewOverlap(si, sj)) {
          // Report new shape overlap
          let e = this.beginContactEvent;
          e.shapeA = si;
          e.shapeB = sj;
          e.bodyA = bi;
          e.bodyB = bj;
          // Reset contact equations
          e.contactEquations.length = 0;
          if (typeof (numContacts) === "number") {
            for (let i = np.contactEquations.length - numContacts; i < np.contactEquations.length; i++) {
              e.contactEquations.push(np.contactEquations[i]);
            }
          }
          this.emit(e);
        }
        // divide the max friction force by the number of contacts
        if (typeof (numContacts) === "number" && numFrictionEquations > 1) { // Why divide by 1?
          for (let i = np.frictionEquations.length - numFrictionEquations; i < np.frictionEquations.length; i++) {
            let f = np.frictionEquations[i];
            f.setSlipForce(f.getSlipForce() / numFrictionEquations);
          }
        }
      }
    }
  }
  /**Add a spring to the simulation
  *
  * @param {Spring} spring
  */
  World.prototype.addSpring(spring) {
    this.springs.push(spring);
    let evt = this.addSpringEvent;
    evt.spring = spring;
    this.emit(evt);
    evt.spring = null;
  }
  /**Remove a spring
  *
  * @param {Spring} spring
  */
  World.prototype.removeSpring(spring) {
    let idx = this.springs.indexOf(spring);
    if (idx !== -1) {
      Utils.splice(this.springs, idx, 1);
    }
  }
  /**Add a body to the simulation
  *
  * @param {Body} body
  *
  * @example
  *     let world = new World(),
  *         body = new Body();
  *     world.addBody(body);
  * @todo What if this is done during step?
  */
  World.prototype.addBody(body) {
    if (this.bodies.indexOf(body) === -1) {
      this.bodies.push(body);
      body.world = this;
      let evt = this.addBodyEvent;
      evt.body = body;
      this.emit(evt);
      evt.body = null;
    }
  }
  /**Remove a body from the simulation. If this method is called during step(), the body removal is scheduled to after the step.
  *
  * @param {Body} body
  */
  World.prototype.removeBody(body) {
    if (this.stepping) {
      this.bodiesToBeRemoved.push(body);
    } else {
      body.world = null;
      let idx = this.bodies.indexOf(body);
      if (idx !== -1) {
        Utils.splice(this.bodies, idx, 1);
        this.removeBodyEvent.body = body;
        body.resetConstraintVelocity();
        this.emit(this.removeBodyEvent);
        this.removeBodyEvent.body = null;
      }
    }
  }
  /**Get a body by its id.
  * @param {number} id
  * @returns {Body} The body, or false if it was not found.
  */
  World.prototype.getBodyById(id) {
    let bodies = this.bodies;
    for (let i = 0; i < bodies.length; i++) {
      let b = bodies[i];
      if (b.id === id) {
        return b;
      }
    }
    return false;
  }
  /**Disable collision between two bodies
  * @param {Body} bodyA
  * @param {Body} bodyB
  */
  World.prototype.disableBodyCollision(bodyA, bodyB) {
    this.disabledBodyCollisionPairs.push(bodyA, bodyB);
  }
  /**Enable collisions between the given two bodies
  * @param {Body} bodyA
  * @param {Body} bodyB
  */
  World.prototype.enableBodyCollision(bodyA, bodyB) {
    let pairs = this.disabledBodyCollisionPairs;
    for (let i = 0; i < pairs.length; i += 2) {
      if ((pairs[i] === bodyA && pairs[i + 1] === bodyB) || (pairs[i + 1] === bodyA && pairs[i] === bodyB)) {
        pairs.splice(i, 2);
        return;
      }
    }
  }
  /**Resets the World, removes all bodies, constraints and springs.
  *
  */
  World.prototype.clear() {
    this.time = 0;
    // Remove all solver equations
    if (this.solver && this.solver.equations.length) {
      this.solver.removeAllEquations();
    }
    // Remove all constraints
    let cs = this.constraints;
    for (let i = cs.length - 1; i >= 0; i--) {
      this.removeConstraint(cs[i]);
    }
    // Remove all bodies
    let bodies = this.bodies;
    for (let i = bodies.length - 1; i >= 0; i--) {
      this.removeBody(bodies[i]);
    }
    // Remove all springs
    let springs = this.springs;
    for (let i = springs.length - 1; i >= 0; i--) {
      this.removeSpring(springs[i]);
    }
    // Remove all contact materials
    let cms = this.contactMaterials;
    for (let i = cms.length - 1; i >= 0; i--) {
      this.removeContactMaterial(cms[i]);
    }
    World.apply(this);
  }
  let hitTest_tmp1 = vec2.create(),
    hitTest_zero = vec2.fromValues(0, 0),
    hitTest_tmp2 = vec2.fromValues(0, 0);
  /**Test if a world point overlaps bodies
  * @param {Array}  worldPoint  Point to use for intersection tests
  * @param {Array}  bodies      A list of objects to check for intersection
  * @param {Number} precision   Used for matching against particles and lines. Adds some margin to these infinitesimal objects.
  * @returns {Array}              Array of bodies that overlap the point
  * @todo Should use an api similar to the raycast function
  * @todo Should probably implement a .containsPoint method for all shapes. Would be more efficient
  * @todo Should use the broadphase
  */
  World.prototype.hitTest(worldPoint, bodies, precision) {
    precision = precision || 0;
    // Create a dummy particle body with a particle shape to test against the bodies
    let pb = new Body({ position: worldPoint }),
      ps = new Particle(),
      px = worldPoint,
      pa = 0,
      x = hitTest_tmp1,
      zero = hitTest_zero,
      tmp = hitTest_tmp2;
    pb.addShape(ps);
    let n = this.narrowphase,
      result = [];
    // Check bodies
    for (let i = 0, N = bodies.length; i !== N; i++) {
      let b = bodies[i];
      for (let j = 0, NS = b.shapes.length; j !== NS; j++) {
        let s = b.shapes[j];
        // Get shape world position + angle
        vec2.rotate(x, s.position, b.angle);
        vec2.add(x, x, b.position);
        let a = s.angle + b.angle;
        if ((s instanceof Circle && n.circleParticle(b, s, x, a, pb, ps, px, pa, true)) ||
          (s instanceof Convex && n.particleConvex(pb, ps, px, pa, b, s, x, a, true)) ||
          (s instanceof Plane && n.particlePlane(pb, ps, px, pa, b, s, x, a, true)) ||
          (s instanceof Capsule && n.particleCapsule(pb, ps, px, pa, b, s, x, a, true)) ||
          (s instanceof Particle && vec2.squaredLength(vec2.sub(tmp, x, worldPoint)) < precision * precision)
        ) {
          result.push(b);
        }
      }
    }
    return result;
  }
  /**Set the stiffness for all equations and contact materials.
  * @param {Number} stiffness
  */
  World.prototype.setGlobalStiffness(stiffness) {
    // Set for all constraints
    let constraints = this.constraints;
    for (let i = 0; i !== constraints.length; i++) {
      let c = constraints[i];
      for (let j = 0; j !== c.equations.length; j++) {
        let eq = c.equations[j];
        eq.stiffness = stiffness;
        eq.needsUpdate = true;
      }
    }
    // Set for all contact materials
    let contactMaterials = this.contactMaterials;
    for (let i = 0; i !== contactMaterials.length; i++) {
      let c = contactMaterials[i];
      c.stiffness = c.frictionStiffness = stiffness;
    }
    // Set for default contact material
    let c = this.defaultContactMaterial;
    c.stiffness = c.frictionStiffness = stiffness;
  }
  /**Set the relaxation for all equations and contact materials.
  * @param {Number} relaxation
  */
  World.prototype.setGlobalRelaxation(relaxation) {
    // Set for all constraints
    for (let i = 0; i !== this.constraints.length; i++) {
      let c = this.constraints[i];
      for (let j = 0; j !== c.equations.length; j++) {
        let eq = c.equations[j];
        eq.relaxation = relaxation;
        eq.needsUpdate = true;
      }
    }
    // Set for all contact materials
    for (let i = 0; i !== this.contactMaterials.length; i++) {
      let c = this.contactMaterials[i];
      c.relaxation = c.frictionRelaxation = relaxation;
    }
    // Set for default contact material
    let c = this.defaultContactMaterial;
    c.relaxation = c.frictionRelaxation = relaxation;
  }
  let tmpAABB = new AABB();
  let tmpArray = [];
  /**Ray cast against all bodies in the world.
  * @param {RaycastResult} result
  * @param {Ray} ray
  * @returns {Boolean} True if any body was hit.
  *
  * @example
  *     let ray = new Ray({
  *         mode: Ray.CLOSEST, // or ANY
  *         from: [0, 0],
  *         to: [10, 0],
  *     });
  *     let result = new RaycastResult();
  *     world.raycast(result, ray);
  *
  *     // Get the hit point
  *     let hitPoint = vec2.create();
  *     result.getHitPoint(hitPoint, ray);
  *     console.log('Hit point: ', hitPoint[0], hitPoint[1], ' at distance ' + result.getHitDistance(ray));
  *
  * @example
  *     let ray = new Ray({
  *         mode: Ray.ALL,
  *         from: [0, 0],
  *         to: [10, 0],
  *         callback: function(result){
  *
  *             // Print some info about the hit
  *             console.log('Hit body and shape: ', result.body, result.shape);
  *
  *             // Get the hit point
  *             let hitPoint = vec2.create();
  *             result.getHitPoint(hitPoint, ray);
  *             console.log('Hit point: ', hitPoint[0], hitPoint[1], ' at distance ' + result.getHitDistance(ray));
  *
  *             // If you are happy with the hits you got this far, you can stop the traversal here:
  *             result.stop();
  *         }
  *     });
  *     let result = new RaycastResult();
  *     world.raycast(result, ray);
  */
  World.prototype.raycast(result, ray) {
    // Get all bodies within the ray AABB
    ray.getAABB(tmpAABB);
    this.broadphase.aabbQuery(this, tmpAABB, tmpArray);
    ray.intersectBodies(result, tmpArray);
    tmpArray.length = 0;
    return result.hasHit();
  }
}, { "../../package.json": 6, "../collision/AABB": 7, "../collision/Broadphase": 8, "../collision/Narrowphase": 10, "../collision/Ray": 11, "../collision/SAPBroadphase": 13, "../constraints/Constraint": 14, "../constraints/DistanceConstraint": 15, "../constraints/GearConstraint": 16, "../constraints/LockConstraint": 17, "../constraints/PrismaticConstraint": 18, "../constraints/RevoluteConstraint": 19, "../events/EventEmitter": 26, "../material/ContactMaterial": 27, "../material/Material": 28, "../math/vec2": 30, "../objects/Body": 31, "../objects/LinearSpring": 32, "../objects/RotationalSpring": 33, "../shapes/Capsule": 38, "../shapes/Circle": 39, "../shapes/Convex": 40, "../shapes/Line": 42, "../shapes/Particle": 43, "../shapes/Plane": 44, "../shapes/Shape": 45, "../solver/GSSolver": 46, "../solver/Solver": 47, "../utils/OverlapKeeper": 52, "../utils/Utils": 57, "./IslandManager": 59 }]
  }, { }, [36])
(36)
});