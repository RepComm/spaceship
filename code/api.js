
export default class API {
  constructor () {
    /**@type {import("../modules/p2.js")}*/
    this.p2;
    /**@type {p2.World} */
    this.pworld;
    /**@type {CanvasRenderingContext2D} */
    this.ctx;

    /**@type {import("./input.js").InputManager} */
    this.input;
  }
}
