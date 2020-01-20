
import { rect, on } from "./aliases.js";
import { Utils, radians, ndist, dist, magnitude, direction } from "./math.js";
import { Entity } from "./entity/entity.js";

class Renderer {
    /**
     * @param {HTMLCanvasElement} canvas renderer element
     * @param {Boolean} renderGrid should the grid render or not
     * @param {Number} gridSpacing coordinate spacing between grid vertices
     */
    constructor(canvas, renderGrid = true, gridSpacing = 1) {
        this.canvas = canvas;
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

    addEntity (entity) {
        this.entities.push (entity);
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
        
        let d, diffX, diffY, mag, dSquared;
        for (let entity of this.entities) {
            entity.step();
            for (let other of this.entities) {
                if (entity === other) continue;
                d = dist(entity.x, entity.y, other.x, other.y);
                dSquared = d*d;
                
                diffX = entity.x - other.x;
                diffY = entity.y - other.y;
                mag = magnitude(diffX, diffY);
                if (mag !== 0) {
                    diffX = diffX / mag;
                    diffY = diffY / mag;
                }
                if (d !== 0) {
                    if (entity.isDynamic) {
                        if (entity.collider && other.collider) {
                            entity.collider.collidesWith(other.collider);
                            if (entity.collider.info.hit) {
                                entity.collider.info.angle -= Math.atan2(entity.velocity.x, entity.velocity.y);
                                entity.velocity.x = (entity.velocity.abs*0.6) * Math.cos(entity.collider.info.angle);
                                entity.velocity.y = (entity.velocity.abs*0.6) * Math.sin(entity.collider.info.angle);
                                
                            } else {
                                entity.addForce(
                                    -(diffX) / dSquared,
                                    -(diffY) / dSquared
                                );
                            }
                        }
                    } else {
                        entity.addForce(
                            -(diffX) / dSquared,
                            -(diffY) / dSquared
                        );
                    }
                }
            }
            this.ctx.save();
            this.ctx.translate(entity.x, entity.y);
            this.ctx.rotate(entity.rotation);
            this.ctx.scale(entity.scale, entity.scale);
            entity.render(this.ctx);
            this.ctx.restore();

            
            d = dist(entity.x, entity.y, this.cursor.localx, this.cursor.localy);
            if (d < 5000) {
                if (d > 5) {
                    this.ctx.beginPath();
                    this.ctx.moveTo(entity.x, entity.y);
                    //this.ctx.lineWidth = (this.gridLineWidth * 4 / this.zoom) / (d/200);
                    this.ctx.lineTo(this.cursor.localx, this.cursor.localy);
                    this.ctx.strokeStyle = "white";
                    this.ctx.stroke();
                } else {
                    this.ctx.fillStyle = "white";
                    this.ctx.save();
                    this.ctx.scale(1/this.zoom, 1/this.zoom);
                    this.ctx.fillText(
                        "X:" + entity.x.toFixed(2) +
                        ", Y:" + entity.y.toFixed(2) +
                        ", V:" + entity.velocity.abs.toFixed(2),
                        this.cursor.localx * this.zoom,
                        this.cursor.localy * this.zoom
                    );
                    this.ctx.restore();
                }
            }
        }

        this.ctx.beginPath();
        this.ctx.ellipse(this.cursor.localx, this.cursor.localy, 1, 1, 0, 0, Math.PI*2);
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
