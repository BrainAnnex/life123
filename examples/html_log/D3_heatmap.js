// A simple D3-based function to draw a heatmap

function draw_heatmap(svg_id, dt)
/*  Expects a div tag
    svg_id: the name of the ID of the SVG element to place the plot into
 */
{
    console.log("Inside draw_heatmap")

    // Labels of row and columns
    const myGroups = ["A", "B", "C"];
    const myVars = ["v1", "v2"];

    // Data for the heatmap
    const MY_DATA = dt;


    // set the dimensions and margins of the graph
    const margins = {top: 30, right: 30, bottom: 30, left: 30};


    const width = 450 - margins.left - margins.right,
          height = 450 - margins.top - margins.bottom;


    // append the svg object to the body of the page
    const svg = d3.select(`#${svg_id}`)
            .append("svg")
                .attr("width", width + margins.left + margins.right)
                .attr("height", height + margins.top + margins.bottom)
            .append("g")
                .attr("transform", `translate(${margins.left},${margins.top})`);


    // Build X scale and axis
    const x = d3.scaleBand()
            .range([ 0, width ])
            .domain(myGroups);

    svg.append("g")
        .attr("transform", `translate(0, ${height})`)
        .call(d3.axisBottom(x));

    // Build Y scale and axis
    const y = d3.scaleBand()
      .range([ height, 0 ])
      .domain(myVars);

    svg.append("g")
      .call(d3.axisLeft(y));

    // Build color scale
    const myColor = d3.scaleLinear()
      .range(["white", "#69b3a2"])
      .domain([1,100]);


    // Transform the DOM element with the SVG
    svg.selectAll()
            .data(MY_DATA, function(d) {return d.group+':'+d.variable;})
            .join("rect")
            .attr("x", function(d) { return x(d.group) })
            .attr("y", function(d) { return y(d.variable) })
            .attr("width", x.bandwidth() )
            .attr("height", y.bandwidth() )
            .style("fill", function(d) { return myColor(d.value)} );
}