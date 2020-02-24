
import { InputManager } from "../input.js";
import { Ship } from "../entity/ship.js";

class ShipController {
  /**
   * @param {InputManager} input 
   * @param {Ship} ship 
   */
  constructor(api, ship) {
    this.api = api;
    this.input = this.api.input;
    this.ship = ship;
    this.rotationalAccelleration = 0.05;
    this.acelleration = 200;
    this.VEC2_BACKWARD = [-this.acelleration, 0];
    this.VEC2_FORWARD = [this.acelleration, 0];
    this.VEC2_CENTER = [0, 0];
  }

  step() {
    if (this.input.isKeyPressed("a")) {
      this.ship.body.angularVelocity = 0;
      this.ship.body.angle -= this.rotationalAccelleration;
    } else if (this.input.isKeyPressed("d")) {
      this.ship.body.angularVelocity = 0;
      this.ship.body.angle += this.rotationalAccelleration;
    }
    if (this.input.isKeyPressed("w")) {
      this.ship.body.applyForceLocal(this.VEC2_FORWARD);
    } else if (this.input.isKeyPressed("s")) {
      //this.ship.body.applyForceLocal(this.VEC2_BACKWARD);
      this.ship.body.velocity[0] *= 0.9;
      this.ship.body.velocity[1] *= 0.9;
    }
  }
}

export { ShipController };
