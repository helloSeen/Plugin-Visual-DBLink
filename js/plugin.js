function FilterDropDown() {

    //Get the value of exact match checkbox
    var checkBox = document.getElementById("myCheck");

    //Check to see if checkbox is clicked 
    if (myCheck.checked == true)
    {
        //Filter table based on the exact coverages
    }

    else
    {
        //Grab the user input min PCT Cov and Q Cov values       
        var PCT = document.getElementById("minPCT").value;
        var Q = document.getElementById("minQ").value;

        var results = document.getElementById("results");

        //Error check the input values
        if (PCT.length == 0 && Q.length == 0)
        {
            alert("Set either a min PCT cov or min Q cov value before clicking the filter button")
            return
        }

        else if (PCT > 100 || PCT < 0)
        {
            alert("Set the minimum PCT coverage value between 0-100");
            return
        }

        else if  (Q > 100 || Q < 0)
        {
            alert("Set the minimum Q coverage value between 0-100");
            return
        }

        //Filter table based on coverages provided

        //Code that adds elements to the dropdown
        //var ele = document.createElement("a");
        //ele.classList = "dropdown-item";
        //ele.href = "#";
        //ele.innerText = "Rahul Varki"; //Replace with whatever we want 
        //document.querySelector(".dropdown-menu").appendChild(ele);
        
    }
}

