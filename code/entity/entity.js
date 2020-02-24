
class Entity {
  /**@param {import("../api.js").default} api*/
  constructor(api) {
    this.api = api;
    /**@type {Path2D} */
    this.path;
    /**@type {p2.Body} */
    this.body;
    /**@type {Boolean} */
    this.isGravitySource = false;
  }

  /**@param {Path2D} path*/
  setPath(path) {
    this.path = path;
  }

  /**
   * @param {CanvasRenderingContext2D} ctx 
   */
  render(ctx=this.api.ctx) {
  }

  get x () {
    return this.body.position[0];
  }

  set x (v) {
    this.body.position[0] = v;
  }

  get y () {
    return this.body.position[1];
  }

  set y (v) {
    this.body.position[1] = v;
  }
}

export { Entity };
