
import {Collider} from "../colliders.js";

class Entity {
    constructor () {
        this.x = 0;
        this.y = 0;
        this.rotation = 0;
        this.velocity = {x:0, y:0, abs:0};
        this.rotationalVelocity = 0;
        /**@type {Path2D} */
        this.path;
        this.mass = 0.1;
        this.isDynamic = true;
        this.friction = 0;
        this.massIsScale = false;
        /**@type {Collider} */
        this.collider;
        this.scale = 1;
    }

    addForce (x, y) {
        this.velocity.x += x;
        this.velocity.y += y;
        this.velocity.abs = Math.sqrt(Math.pow(this.velocity.x, 2) + Math.pow(this.velocity.y, 2));
    }

    mulForce (val) {
        this.velocity.x *= val;
        this.velocity.y *= val;
        this.velocity.abs = Math.sqrt(Math.pow(this.velocity.x, 2) + Math.pow(this.velocity.y, 2));
    }

    step () {
        if (this.isDynamic) {
            this.rotation += this.rotationalVelocity / this.mass;
            this.x += this.velocity.x / this.mass;
            this.y += this.velocity.y / this.mass;
        }
        this.collider.x = this.x;
        this.collider.y = this.y;
    }

    setPath (path) {
        this.path = path;
    }

    render (ctx) {
    }
}

export {Entity};
