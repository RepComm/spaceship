
import {dist, direction} from "./math.js";

class Collider {
    constructor () {
        this.type = "collider";
        this.x = 0;
        this.y = 0;
        this.info = {
            distance:Infinity, //Distance of centers
            hit:false, //Hit or not
            point:{x:0, y:0} //Hit point
        };
    }

    /**Check collision info of this one and another shape
     * @param {Collider} other 
     */
    collidesWith (other) {

    }
}

class Circle extends Collider {
    constructor () {
        super();
        this.type = "circle";
        this.radius = 1;
    }

    /**Check collision info of this one and another shape
     * @param {Collider} other 
     */
    collidesWith (other) {
        switch(other.type) {
            case Circle.type:
                this.info.distance = dist(this.x, this.y, other.x, other.y);
                this.info.hit = this.info.distance < this.radius + other.radius;
                this.info.point = direction(this.x, this.y, other.x, other.y);
                this.info.point.x *= this.radius;
                this.info.point.y *= this.radius;
                this.info.angle = Math.atan2(this.y - other.y, this.x - other.x);
                break;
            default:
                break;
        }
        //return this.info;
    }
}
Circle.type = "circle";

export {Circle, Collider};
