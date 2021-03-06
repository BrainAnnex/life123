// A simple D3-based function to draw a bar chart

function draw_barchart(svg_id, dt)
/*  Expects a svg tag
    svg_id: the name of the ID of the SVG element to place the plot into
 */
{
    console.log("Inside draw_barchart")

    const container = d3.select(`#${svg_id}`)
        .classed('container', true);


    const bars = container
        .selectAll('.bar')
        .data(dt)
        .enter()
        .append('rect')
        .classed('bar', true)
        .attr('width', '54')
        .attr('height', dt => (dt.value * 15))
        .attr('x', dt => 8 + (dt.index * 60))
        .attr('y', dt => 200 - (dt.value * 15));
}