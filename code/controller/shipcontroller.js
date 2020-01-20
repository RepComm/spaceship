
import { InputManager } from "../input.js";
import { Ship } from "../entity/ship.js";
import { magnitude } from "../math.js";

class ShipController {
    /**
     * @param {InputManager} input 
     * @param {Ship} ship 
     */
    constructor(input, ship) {
        this.input = input;
        this.ship = ship;
        this.rotationalAccelleration = 0.01;
        this.acelleration = 0.1;
    }

    step() {
        if (this.input.isKeyPressed("a")) {
            //this.ship.rotationalVelocity -= this.rotationalAccelleration;
            this.ship.rotation -= 0.05;
        } else if (this.input.isKeyPressed("d")) {
            //this.ship.rotationalVelocity += this.rotationalAccelleration;
            this.ship.rotation += 0.05;
        }
        if (this.input.isKeyPressed("w")) {
            this.ship.addForce(
                Math.cos(this.ship.rotation) * this.acelleration,
                Math.sin(this.ship.rotation) * this.acelleration
            );
        } else if (this.input.isKeyPressed("s")) {
            this.ship.mulForce(0.95);
        }
    }
}

export { ShipController };
