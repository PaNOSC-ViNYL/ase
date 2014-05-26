function open_row(x, n) {
    r = x.rowIndex;
    request = XMLHttpRequest();
    request.open('GET', '/open_row/' + n, true);
    request.onload = function() {
        data = request.responseText;
        var table = document.getElementById('rows');
        if (data) {
            var row = table.insertRow(r + 1);
            var cell = row.insertCell(0);
            cell.colSpan = 100;
            cell.innerHTML = data;
        } else {
            table.deleteRow(r + 1);
        }
    }
    request.send();
}
