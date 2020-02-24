
import { rect, on } from "./aliases.js";
import { Utils, dist, magnitude, direction } from "./math.js";
import { Entity } from "./entity/entity.js";
import API from "./api.js";

class Renderer {
  /**@param {API} api
   * @param {Boolean} renderGrid should the grid render or not
   * @param {Number} gridSpacing coordinate spacing between grid vertices
   */
  constructor(api, renderGrid = true, gridSpacing = 1) {
    this.api = api;
    this.canvas = this.api.canvas;
    this.ctx = canvas.getContext("2d");

    this.renderGrid = renderGrid;
    this.gridSpacing = gridSpacing;
    this.gridColor = "white";
    this.gridLineWidth = 0.1;

    this.centerX = 0;
    this.centerY = 0;
    this.left = 0;
    this.right = 0;
    this.top = 0;
    this.bottom = 0;

    this.zoom = 100;

    this.drawRect;
    this.prevDrawRect = { width: 0, height: 0 };

    this.cursor = { x: 0, y: 0, localx: 0, localy: 0, };

    /**@type {Array<Entity>} */
    this.entities = new Array();

    this.renderRequestCallback = () => {
      this.render();
    };

    on(window, "resize", () => {
      console.log("Resize");
    });

    requestAnimationFrame(this.renderRequestCallback);
  }

  addEntity(entity) {
    this.entities.push(entity);
  }

  /** Set the page-space (layerX/layerY) cursor position
   * Triggers a render
   * @param {Integer} x 
   * @param {Integer} y
   */
  setCursor(x, y) {
    this.cursor.x = x;
    this.cursor.y = y;
    //this.needsRender = true;
  }

  /** Set the origin the viewer is viewing from
   * Triggers a render
   * @param {Integer} x
   * @param {Integer} y 
   */
  setCenter(x, y) {
    this.centerX = x;
    this.centerY = y;
  }

  /** Move the origin the viewer is viewing from by some amounts
   * @param {Integer} xa move amount x
   * @param {Integer} ya move amount y
   */
  moveCenter(xa, ya) {
    this.setCenter(this.centerX + xa, this.centerY + ya);
  }

  /** Trigger a render
   * NOTE: May not actually render when <Renderer>.needsRender not true
   */
  render() {
    requestAnimationFrame(this.renderRequestCallback);

    this.api.pworld.step(1 / 30, 1 / 30, 10);

    this.drawRect = rect(this.canvas);
    if (this.prevDrawRect.width !== this.drawRect.width ||
      this.prevDrawRect.height !== this.drawRect.height) {
      this.prevDrawRect.width = this.drawRect.width;
      this.prevDrawRect.height = this.drawRect.height;
      this.canvas.width = this.drawRect.width;
      this.canvas.height = this.drawRect.height;
    }

    this.left = (-this.drawRect.width / 2) / this.zoom;
    this.right = (this.drawRect.width / 2) / this.zoom;
    this.top = (-this.drawRect.height / 2) / this.zoom;
    this.bottom = (this.drawRect.height / 2) / this.zoom;

    this.ctx.save();
    this.ctx.clearRect(0, 0, this.drawRect.width, this.drawRect.height);
    this.ctx.translate(this.drawRect.width / 2, this.drawRect.height / 2);

    this.ctx.scale(this.zoom, this.zoom);
    this.ctx.lineWidth = this.gridLineWidth / this.zoom;
    this.ctx.translate(-this.centerX, -this.centerY);

    if (this.renderGrid) {
      this.ctx.beginPath();
      let xOffset = Utils.roundTo(this.left + this.centerX, this.gridSpacing);
      for (let x = xOffset; x < this.right + this.centerX; x += this.gridSpacing) {
        this.ctx.moveTo(x, this.top + this.centerY);
        this.ctx.lineTo(x, this.bottom + this.centerY);
      }
      let yOffset = Utils.roundTo(this.top + this.centerY, this.gridSpacing);
      for (let y = yOffset; y < this.bottom + this.centerY; y += this.gridSpacing) {
        this.ctx.moveTo(this.left + this.centerX, y);
        this.ctx.lineTo(this.right + this.centerX, y);
      }
      this.ctx.strokeStyle = this.gridColor;
      this.ctx.stroke();
    }

    this.ctx.beginPath();

    this.cursor.localx = ((this.cursor.x - this.drawRect.width / 2) / this.zoom) + this.centerX;
    this.cursor.localy = ((this.cursor.y - this.drawRect.height / 2) / this.zoom) + this.centerY;

    let d, dir;
    for (let entity of this.entities) {
      for (let other of this.entities) {
        if (other !== entity) {
          if (other.isGravitySource) {
            d = dist(entity.x, entity.y, other.x, other.y) * 2;
            dir = direction(entity.x, entity.y, other.x, other.y);
            entity.body.applyForce([
              (Math.cos(dir) * entity.body.mass) / d, (Math.sin(dir) * entity.body.mass) / d
            ]);
          }
        }
      }
      this.ctx.save();
      this.ctx.translate(entity.body.position[0], entity.body.position[1]);
      this.ctx.rotate(entity.body.angle);
      this.ctx.scale(entity.scale, entity.scale);
      entity.render(this.ctx);
      this.ctx.restore();
    }

    this.ctx.beginPath();
    this.ctx.ellipse(this.cursor.localx, this.cursor.localy, 1, 1, 0, 0, Math.PI * 2);
    this.ctx.closePath();
    this.ctx.strokeStyle = "white";
    this.ctx.stroke();

    this.ctx.restore();
  }

  /** Add some zoom to your room..
   * Triggers a render
   * @param {Number} za zoom amount
   */
  addZoom(za) {
    this.zoom -= za;
    if (this.zoom < 8) {
      this.zoom = 8;
    } else if (this.zoom > 400) {
      this.zoom = 400;
    }
  }
}

export { Renderer };
